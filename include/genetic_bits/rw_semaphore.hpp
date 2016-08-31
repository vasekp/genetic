namespace gen {

/* Internal functions only to be used from other sources in this directory. */
namespace internal {

  /* Writers preference */
  struct rw_semaphore {
    protected:
    std::mutex m{};  // protects both the members of this and the resource
    unsigned r{0}, w{0};  // active readers, waiting writers
    std::condition_variable rq{}, wq{};  // read queue, write queue
    size_t mod_cnt{0};  // modification counter

    friend class rw_lock;

    public:
    rw_semaphore() { }

    rw_semaphore(const rw_semaphore&) = delete;
    rw_semaphore(rw_semaphore&&) = delete;
    rw_semaphore& operator= (const rw_semaphore&) = delete;
    rw_semaphore& operator= (rw_semaphore&&) = delete;

    size_t get_mod_cnt() { return mod_cnt; }
  };


  class rw_lock {
    /* Invariants:
     *  mutex free ⇒ semaphore members and resource constant
     *  pers_lock free ⇒ no active writes (resource constant)
     *  w = 0 ⇒ rq notified
     *  r = 0 ⇒ wq notified
     *  write active ⇒ r = 0, w > 0, pers_lock */
    rw_semaphore& s;
    std::unique_lock<std::mutex> pers_lock{};
    bool write{false};

    protected:
    rw_lock(rw_semaphore& _s, bool _w): s(_s), write(_w) {
      if(write) {
        // Logic: write creates a persistent lock, preventing other writes from
        // starting. They can still execute ++s.w to let themselves known but
        // this enqueues them at the end of wq.wait. Write can only start when
        // everyone is finished reading.
        pers_lock = std::unique_lock<std::mutex>(s.m);
        ++s.w;
        s.wq.wait(pers_lock, [&] { return s.r == 0; });
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
        ++s.mod_cnt;
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
    bool upgrade() {
      if(write)
        return false;
      // Atomically combined action of -read and +write
      pers_lock = std::unique_lock<std::mutex>(s.m);
      ++s.w;
      --s.r;
      write = true;
      s.wq.wait(pers_lock, [&] { return s.r == 0; });
      // It makes no sense for another waiting write to wake up now, that
      // would defeat the purpose of the upgrade. But we need to notify the
      // other threads that r might be zero as no one else would do that.
      // Keeping pers_lock allows them to exit wait() but stay stuck at the
      // mutex lock.
      s.wq.notify_all();
      return true;
    }

    void downgrade() {
      if(!write)
        return;
      // Atomically combined action of -write and +read.
      --s.w;
      ++s.mod_cnt;
      ++s.r;
      write = false;
      pers_lock.unlock();
      // Other waiting reads can proceed now, the write is done
      s.rq.notify_all();
    }
  };


  class read_lock: public rw_lock {
    public:
    read_lock(rw_semaphore& s): rw_lock(s, false) { }
  };

  class write_lock: public rw_lock {
    public:
    write_lock(rw_semaphore& s): rw_lock(s, true) { }
  };

  class upgrade_lock {
    rw_lock& l;
    bool u{false};

    public:
    upgrade_lock(rw_lock& _l): l(_l) {
      u = l.upgrade();
    }

    ~upgrade_lock() {
      if(u)
        l.downgrade();
    }
  };

} // namespace internal

} // namespace gen
