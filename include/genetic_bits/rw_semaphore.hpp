namespace gen {

/* Internal functions only to be used from other sources in this directory. */
namespace internal {

  /*** Writers' preference ***/
  /* Invariants:
   *  mutex free ⇒ semaphore members and resource constant
   *  pers_lock free ⇒ no active writes (resource constant)
   *  w = 0 ⇔ no writers waiting
   *  r = 0 ⇔ no readers active
   *  w → 0 ⇒ rq notified
   *  r → 0 ⇒ wq notified
   *  cnt = cache_cnt ⇒ res = cache_res
   *  write active ⇒ r = 0, w > 0, pers_lock, cnt > cache_cnt
   * Corollary:
   *  see rw_lock::upgrade()! */
  class rw_semaphore {

  protected:

    std::mutex m{};  // protects both the members of this and the resource
    unsigned r{0}, w{0};  // active readers, waiting writers
    std::condition_variable rq{}, wq{};  // read queue, write queue
    size_t mod_cnt{0};  // modification counter

    friend class rw_lock;

  public:

    rw_semaphore() { }

    /* Copy: creates a new semaphore but keeps mod_cnt
     * (so that if (x_,s_) are copies of (x,s) where x.last_mod_cnt =
     * s.get_mod_cnt() then x_.last_mod_cnt = s_.get_mod_cnt() */
    rw_semaphore(const rw_semaphore& s): mod_cnt(s.mod_cnt) { }

    /* Move: everything can be reused */
    rw_semaphore(rw_semaphore&&) = default;

    rw_semaphore& operator= (const rw_semaphore& s) {
      mod_cnt = s.mod_cnt;
      return *this;
    }

    rw_semaphore& operator= (rw_semaphore&&) = default;

    ~rw_semaphore() noexcept(false) {
      if(r != 0 || w != 0)
        throw std::logic_error("Semaphore destroyed in use!");
    }

    size_t get_mod_cnt() { return mod_cnt; }

  }; // class rw_semaphore


  class rw_lock {

    rw_semaphore& s;
    std::unique_lock<std::mutex> pers_lock{};
    bool write{false};

  protected:

    rw_lock(rw_semaphore& s_, bool write_, bool inc): s(s_), write(write_) {
      if(write) {
        // Logic: write creates a persistent lock, preventing other writes from
        // starting. They can still execute ++s.w to let themselves known but
        // this enqueues them at the end of wq.wait. Write can only start when
        // everyone is finished reading.
        pers_lock = std::unique_lock<std::mutex>(s.m);
        ++s.w;
        s.wq.wait(pers_lock, [&] { return s.r == 0; });
        // Increment the modification counter conditionally. Sometimes a
        // writer lock is acquired and later found that no modification is
        // necessary, e.g., if another thread did it in the meantime. Another
        // situation where this can be useful is when an update of the
        // resource does not count as a change of the observable state.
        if(inc)
          bump();
      } else {
        // Logic: read creates a temporary lock for the purposes of rq.wait and
        // s.r access. Read can only start when no one is writing or waiting to
        // write.
        std::unique_lock<std::mutex> temp_lock(s.m);
        s.rq.wait(temp_lock, [&] { return s.w == 0; });
        ++s.r;
      }
    }

    ~rw_lock() {
      if(write) {
        --s.w;
        pers_lock.unlock();
        s.rq.notify_all();
      } else {
        {
          std::unique_lock<std::mutex> temp_lock(s.m);
          --s.r;
        }
        s.wq.notify_all();
      }
    }

    rw_lock(const rw_lock&) = delete;
    rw_lock(rw_lock&&) = delete;
    rw_lock& operator= (const rw_lock&) = delete;
    rw_lock& operator= (rw_lock&&) = delete;

  public:

    /* Upgrade a read lock to a write lock, giving this thread priority over
     * waiting write-only locks.
     *
     * Important: don't rely on any information about the resource obtained
     * prior to upgrade()! If more readers simultaneously want to upgrade then
     * they will have to wait for each other, and may get the lock in any
     * order, changing the state visible to the other ones.
     *
     * Returns whether an upgrade happened: the returned value is false if
     * this was already a writer lock. This is useful for a conditional
     * downgrade() later when given a lock of unknown type. */
    bool upgrade(bool inc = true) {
      if(write)
        return false;
      // Atomically combined action of -read and +write
      pers_lock = std::unique_lock<std::mutex>(s.m);
      ++s.w;
      --s.r;
      write = true;
      s.wq.wait(pers_lock, [&] { return s.r == 0; });
      if(inc)
        bump();
      // It makes no sense for another waiting write to wake up now, that
      // would defeat the purpose of the upgrade. But we need to notify the
      // other threads that r might be zero as no one else would do that.
      // Keeping pers_lock allows them to exit wait() but stay stuck at the
      // mutex lock.
      s.wq.notify_all();
      return true;
    }

    /* Upgrade if a condition is met. If another thread (with an upgraded read
     * lock) breaks the condition in the time of obtaining the write lock,
     * it's released again.
     *
     * Returns whether the upgrade was successful: when the return value is
     * true, we (1) have a write lock and (2) are sure that cond() is true.
     * Throws an exception if this was already a write lock, because it could
     * not be distinguished whether the upgrade was not successful due to the
     * condition or due to that. */
    bool upgrade_if(std::function<bool()> cond, bool inc = true) {
      if(write)
        throw std::logic_error("upgrade_if called on a write lock");
      if(!cond())
        return false;
      upgrade(false);
      if(!cond()) {
        downgrade();
        return false;
      } else {
        if(inc)
          bump();
        return true;
      }
    }

    /* Downgrade a write lock to a read lock, allowing the read to finish
     * before another writer obtains the lock. */
    void downgrade() {
      if(!write)
        return;
      // Atomically combined action of -write and +read.
      --s.w;
      ++s.r;
      write = false;
      pers_lock.unlock();
      // Other waiting reads can proceed now, the write is done
      s.rq.notify_all();
    }

    /* Manually increment the modification counter. This is used in patterns
     * where we don't know for sure whether an actual modification is going to
     * happen at the time of acquiring a writer's lock or upgrading a reader's
     * one. With default parameters to the respective constructors this
     * happens automatically.
     *
     * Can be called multiple times within a write_lock. However, it must be
     * called at least once before any operations affecting the observable
     * state. */
    void bump() {
      if(!write)
        throw std::logic_error("bump() called in a reader lock");
      else
        ++s.mod_cnt;
    }

  }; // class rw_lock


  class read_lock: public rw_lock {

  public:
    read_lock(rw_semaphore& s): rw_lock(s, false, false) { }

  }; // class read_lock


  class write_lock: public rw_lock {

  public:
    write_lock(rw_semaphore& s, bool inc = true): rw_lock(s, true, inc) { }

  }; // class write_lock


  class upgrade_lock {

    rw_lock& l;
    bool u{false};

  public:

    upgrade_lock(rw_lock& l_, bool inc = true): l(l_) {
      u = l.upgrade(inc);
    }

    ~upgrade_lock() {
      if(u)
        l.downgrade();
    }

  }; // class upgrade_lock

} // namespace internal

} // namespace gen
