#ifndef LIBSPANNER_TIMER_H
#define LIBSPANNER_TIMER_H

#include <string>
#include <CGAL/Real_timer.h>

namespace spanner {

    class Timer {
    public:
        explicit Timer(std::string delimiter = ",") : m_delimiter(std::move(delimiter)) {
            m_clock.start();
        }

        ~Timer() {
            if( m_clock.is_running() ) {
                std::cout << stop() << m_delimiter;
            }
        }

        double stop() {
            m_clock.stop();
            return m_clock.time();
        }

    private:
        std::string m_delimiter;
        CGAL::Real_timer m_clock;
    };

} // namespace spanner


#endif // LIBSPANNER_TIMER_H


