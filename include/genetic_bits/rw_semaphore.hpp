namespace gen {

/* Internal functions only to be used from other sources in this directory. */
namespace internal {

  struct rw_semaphore {
    protected:
    unsigned r{0};
    unsigned w{0};
    std::mutex m{};
    std::condition_variable cv{};

    friend class rw_lock;

    public:
    rw_semaphore() { }

    rw_semaphore(const rw_semaphore&) = delete;
    rw_semaphore(rw_semaphore&&) = delete;
    rw_semaphore& operator= (const rw_semaphore&) = delete;
    rw_semaphore& operator= (rw_semaphore&&) = delete;
  };


  class rw_lock {
    rw_semaphore& s;
    std::unique_lock<std::mutex> pers_lock{};
    bool write{false};

    protected:
    rw_lock(rw_semaphore& _s, bool _w): s(_s), write(_w) {
      if(write) {
        // Logic: write creates a persistent lock, preventing other writes from
        // starting. They can still execute ++s.w to let themselves known but
        // this enqueues them at the end of cv.wait. Write can only start when
        // everyone is finished reading.
        pers_lock = std::unique_lock<std::mutex>(s.m);
        ++s.w;
        s.cv.wait(pers_lock, [&] { return s.r == 0; });
      } else {
        // Logic: read creates a temporary lock for the purposes of cv.wait and
        // s.r access. Read can only start when no one is writing or waiting to
        // write.
        std::unique_lock<std::mutex> temp_lock(s.m);
        s.cv.wait(temp_lock, [&] { return s.w == 0; });
        ++s.r;
      }
    }

    ~rw_lock() {
      if(write) {
        --s.w;
        pers_lock.unlock();
      } else {
        std::unique_lock<std::mutex> temp_lock(s.m);
        --s.r;
      }
      // Notify whomever may be waiting for r=0 or w=0
      s.cv.notify_all();
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
      // No cv.notify_all: it makes no sense for a waiting write to wake up now;
      // a read is technically still underway. Moreover, we're a writing thread
      // with work to do.
      pers_lock = std::unique_lock<std::mutex>(s.m);
      ++s.w;
      --s.r;
      write = true;
      s.cv.wait(pers_lock, [&] { return s.r == 0; });
      return true;
    }

    void downgrade() {
      if(!write)
        return;
      // Atomically combined action of -write and +read.
      // cv.notify_all: other waiting reads can proceed now, the write is done
      --s.w;
      ++s.r;
      write = false;
      pers_lock.unlock();
      s.cv.notify_all();
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
