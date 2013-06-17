#ifndef MEAMZ_TIMER_H
#define MEAMZ_TIMER_H

namespace MEAMZ_NS
{

class Timer
{
public:
  Timer();
  virtual ~Timer();

  void start();       // Start timer
  void stop();        // Stop timer
  double elapsed();   // Output elapsed time since start
  double get_time();  // Output the value of the time_ variable

protected:

private:
  double time_;

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TIMER_H
