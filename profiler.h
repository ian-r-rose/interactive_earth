/************************************
  Lightweight phase timers for the
  solver and renderer. Timers are
  scoped, keyed by an enum so the hot
  path never touches strings, and can
  be compiled out with IE_NO_PROFILE.
************************************/
#ifndef PROFILER_H
#define PROFILER_H

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif

enum class Phase
{
  stokes_rhs,
  stokes_fft,
  stokes_tridiag,
  stokes_velocity,
  heat,
  advect,
  diffuse_fft,
  diffuse_tridiag,
  diffuse_misc,
  waves,
  tpw,
  step,  //One full solver timestep
  render_prep,  //CPU-side preparation of vertex data
  render_gl,  //GL calls (CPU submission time only; GL is asynchronous)
  frame,  //Wall time between frames
  count
};

inline const char *phase_name(Phase p)
{
  static const char *names[] = {
    "stokes_rhs", "stokes_fft", "stokes_tridiag", "stokes_velocity",
    "heat", "advect", "diffuse_fft", "diffuse_tridiag", "diffuse_misc",
    "waves", "tpw", "step", "render_prep", "render_gl", "frame" };
  return names[static_cast<int>(p)];
}

//Wall clock in milliseconds
inline double now_ms()
{
#ifdef __EMSCRIPTEN__
  return emscripten_get_now();
#else
  using namespace std::chrono;
  return duration<double, std::milli>(steady_clock::now().time_since_epoch()).count();
#endif
}

class Profiler
{
  public:
    //Keep the most recent `window` samples per phase, or all of them if zero.
    explicit Profiler(size_t window = 120)
      : window(window), samples(n_phases), next(n_phases, 0), pending(n_phases, -1.) {}

    void set_window(size_t w) { window = w; reset(); }

    void reset()
    {
      for (auto &s : samples) s.clear();
      std::fill(next.begin(), next.end(), 0);
      std::fill(pending.begin(), pending.end(), -1.);
    }

    //Accumulate time for a phase in the current step or frame. A phase
    //may be timed in several pieces; commit() records their sum.
    void accumulate(Phase p, double ms)
    {
      double &acc = pending[static_cast<int>(p)];
      acc = (acc < 0. ? ms : acc + ms);
    }

    //Record one sample for every phase accumulated since the last commit.
    void commit()
    {
      for (int i = 0; i < n_phases; ++i)
        if (pending[i] >= 0.)
        {
          record(static_cast<Phase>(i), pending[i]);
          pending[i] = -1.;
        }
    }

    void record(Phase p, double ms)
    {
      const int i = static_cast<int>(p);
      std::vector<double> &s = samples[i];
      if (window == 0 || s.size() < window)
        s.push_back(ms);
      else
      {
        s[next[i]] = ms;
        next[i] = (next[i] + 1) % window;
      }
    }

    bool has(Phase p) const { return !samples[static_cast<int>(p)].empty(); }

    double mean(Phase p) const
    {
      const std::vector<double> &s = samples[static_cast<int>(p)];
      double sum = 0.;
      for (double x : s) sum += x;
      return s.empty() ? 0. : sum / s.size();
    }

    double max(Phase p) const
    {
      const std::vector<double> &s = samples[static_cast<int>(p)];
      return s.empty() ? 0. : *std::max_element(s.begin(), s.end());
    }

    //Percentile in [0,1], nearest rank
    double percentile(Phase p, double q) const
    {
      std::vector<double> s = samples[static_cast<int>(p)];
      if (s.empty()) return 0.;
      std::sort(s.begin(), s.end());
      return s[std::min(s.size() - 1, size_t(q * (s.size() - 1) + 0.5))];
    }

    //JSON object of {phase: {mean, median, p95, max}} for phases with samples
    std::string json() const
    {
      std::string out = "{";
      char buf[160];
      bool first = true;
      for (int i = 0; i < n_phases; ++i)
      {
        const Phase p = static_cast<Phase>(i);
        if (!has(p)) continue;
        std::snprintf(buf, sizeof(buf), "%s\"%s\": {\"mean\": %.4f, \"median\": %.4f, \"p95\": %.4f, \"max\": %.4f}",
                      first ? "" : ", ", phase_name(p), mean(p), percentile(p, 0.5),
                      percentile(p, 0.95), max(p));
        out += buf;
        first = false;
      }
      return out + "}";
    }

  private:
    static const int n_phases = static_cast<int>(Phase::count);
    size_t window;
    std::vector<std::vector<double>> samples;
    std::vector<size_t> next;
    std::vector<double> pending;  //Negative means not timed since last commit
};

//The global profiler used by the solver and renderer
inline Profiler &profiler()
{
  static Profiler instance;
  return instance;
}

class ScopedTimer
{
  public:
    explicit ScopedTimer(Phase p) : phase(p), start(now_ms()) {}
    ~ScopedTimer() { profiler().accumulate(phase, now_ms() - start); }
  private:
    Phase phase;
    double start;
};

#define IE_CONCAT_INNER(a, b) a##b
#define IE_CONCAT(a, b) IE_CONCAT_INNER(a, b)
#ifdef IE_NO_PROFILE
#define IE_PROFILE(phase)
#else
//Time the enclosing scope, e.g. IE_PROFILE(advect);
#define IE_PROFILE(phase) ScopedTimer IE_CONCAT(ie_timer_, __LINE__)(Phase::phase)
#endif

#endif
