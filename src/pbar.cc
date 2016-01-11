/**
 * One-line refreshing progress bar inspired by my usual python pbar.
 * 
 *   description [#---------] k/n  10% [time: 01:23:45s, eta: 01d 00:12:34, 2.7 iters/sec]
 *
 * @author M. Fouesneau
 *  
 *  compile with
 *
 *      g++ -O3 --std=c++11 pbar.cc -o pbar
 *
 *  .. code:: c
 *
 *
 *		size_t n = 10000;
 *		PBar pb(n);
 *      pb.set_description("short loop");
 *
 *      std::cout << "test short loop" << std::endl;
 *      // neat to also increment the pb during the for declaration!!
 *		for(size_t i=0; i <= n; ++i, ++pb) {
 *			usleep(200);
 *		}
 */
#include "pbar.hh"


/** make a human readable tie interval decomposed into days, hours, minutes, and
 * seconds.
 *
 * @param t     interval duration
 * @return txt  string representing the interval (<days>d <hrs>:<min>:<sec>)
 */
std::string PBar::format_interval(duration t){
    typedef std::chrono::duration<size_t, std::ratio<86400>> days;
    using std::chrono::duration_cast;
    using std::chrono::hours;
    using std::chrono::minutes;
    using std::chrono::seconds;

    //percentage
    int d=0, h=0, m=0, s=0;
    
    if ( t >= days(1) ) {
        auto numDays = duration_cast<days>(t);
        t -= numDays;
        d = numDays.count();
    }
    
    if ( t >= hours(1) ) {
        auto numHours = duration_cast<hours>(t);
        t -= numHours;
        h = numHours.count();
    }

    if ( t >= minutes(1) ) { 
        auto numMins = duration_cast<minutes>(t);
        t -= numMins;
        m = numMins.count();
    }
    
    if ( t >= seconds(1) ) {
        s = duration_cast<seconds>(t).count();
    }

    char interval_str[9] = "00:00:00";  // 00:00:00
    sprintf(interval_str, "%02d:%02d:%02d", h, m, s);

    std::string out_str = "";
    if (d > 0){
        out_str.append(std::to_string(d));
        out_str.append("d ");
    }
    out_str.append(interval_str);
    
    return out_str;
}

/**
 * make a progress string.
 *
 * description  [----------] k/n  10% [time: 00:00:00, eta: 00:00:00, 2.7 iters/sec]
 *
 * @param  n     number of finished iterations 
 * @param  total total number of iterations 
 * @param  dt    number of seconds passed since start
 * @return txt   string representing the meter
 * 
 */
std::string PBar::build_str_meter(std::string desc, size_t n, 
        size_t total, duration dt)
{
    using std::chrono::duration_cast;
    using std::chrono::seconds;

    //elapsed time string
    std::string elapsed_str = this->format_interval(dt);
    
    //calculate rate: iter / second
    float rate = static_cast<float>(n) / duration_cast<seconds>(dt).count();
    char rate_str[10];
    sprintf(rate_str, "%5.2f", rate);

    //calculate the bar representation
    float frac = static_cast<float>(n) / static_cast<float>(total);
    size_t bar_length = static_cast<size_t>(frac * this->width);
    std::string bar_str = "|";
    bar_str.append(bar_length, '#');
    bar_str.append(this->width - bar_length, '-');
    bar_str.append("|");

    //percentage
    char percent_str[6];
    sprintf(percent_str, "%3d%%", (int)(100 * frac));

    //ETA, what's left
    std::string eta_str = "?";
    if (n > 0){
        duration esecs = duration_cast<seconds>(dt);
        duration eta = duration_cast<duration>(esecs / n * (total - n));
        eta_str = this->format_interval(eta);
    }

    std::string txt = "";
    txt.append(desc);
    txt.append(" ");

    if (total <= 0){  // total number unknown, no bar
        txt.append(std::to_string(n));
    } else {
        txt.append(bar_str);
        txt.append(" ");
        txt.append(std::to_string(n));
        txt.append("/");
        txt.append(std::to_string(total));
        txt.append(" ");
        txt.append(percent_str);
    }
    // rest of the info
    txt.append(" [time: ");
    txt.append(elapsed_str);
    txt.append(", ");
    if (total > 0){  // add eta
        txt.append("eta: ");
        txt.append(eta_str);
        txt.append(", ");
    }
    txt.append(rate_str);
    txt.append(" ");
    txt.append(this->units);
    txt.append("/sec");
    txt.append("]");

    // pad with empty spaces to clean the line and drop spaces from size
    size_t curlen = txt.length();
    if (txt.length() < this->current_length){
        txt.append(this->current_length - txt.length(), ' ');
    }
    this->current_length = curlen;

    return txt;
}


/** print a status on the last file line and clean the rest of the line.
 *
 * @param s  message to write
 */
void PBar::print_status(std::string s){
    std::cout << '\r';
    std::cout << s;
    std::cout.flush();
}


/** 
 * reset the progress bar.
 *
 * @param _n maximum iteration state
 * */
void PBar::reset( size_t _n=0 ) 
{ 
    this->n = _n;
    this->cur = 0;
    this->desc=""; 
    this->current_length = 0;
    this->units = "iters";
}


/** Start the timer */
void PBar::start() { 
    this->startTime = system_clock::now();
    this->lastCheck = startTime;
}

/** stop the timer */
void PBar::finish() { 
    std::chrono::system_clock::time_point endTime = system_clock::now();
    auto elapsed = endTime - this->startTime;

    // update only at most every 1 second or at the end
    std::string meter_str = this->build_str_meter(this->desc, 
            this->cur,
            this->n, elapsed);
    this->print_status(meter_str);
    this->lastCheck = endTime;
    std::cout << std::endl;
}


/** Set PBar string description.
 *
 * @param desc description string 
 */
void PBar::set_description(std::string desc){
    this->desc = desc;
}


/** Set PBar string units.
 *
 * @param units units description string 
 */
void PBar::set_units(std::string units){
    this->units = units;
}

/** Add progress using operator ++ */
void PBar::operator++() {
    if ((this->cur >= this->n) and not(this -> n <= 0))
        // progress is over already
        return;
    if (this->cur == 0){ // start is not already started
        this->start();
    }

    // update content
    ++this->cur;		
    std::chrono::system_clock::time_point endTime = system_clock::now();

    auto elapsed = endTime - this->startTime;

    // update only at most every 1 second or at the end
    if ( (this->cur <= 1) or 
            (endTime - lastCheck) >= std::chrono::seconds(1) or 
            ((this->cur >= this->n) and (this->n > 0)) ) {
        std::string meter_str = this->build_str_meter(this->desc, 
                this->cur,
                this->n, elapsed);
        this->print_status(meter_str);
        this->lastCheck = endTime;
        if ((this->cur >= this->n) and not(this -> n <= 0)){
            std::cout << std::endl;
        }
    }
}


/** overload +=  to add specific progression value.
 *
 * @param d current iteration
 */
void PBar::operator+=( const size_t d ) {
    // update values
    if (this->cur >= this->n) 
        return;
    if (this->cur == 0){
        this->start();
    }
    this->cur += d;
    if (this->cur > this->n)
        this->cur = this->n;

    std::chrono::system_clock::time_point endTime = system_clock::now();

    auto elapsed = endTime - this->startTime;

    // update only at most every 1 second or at the end
    if ( (endTime - lastCheck) >= std::chrono::seconds(1) or (cur == n) ) {

        std::string meter_str = this->build_str_meter(this->desc, 
                this->cur,
                this->n, elapsed);
        this->print_status(meter_str);
        this->lastCheck = endTime;
    }
}


/*-----------------------------------------------------------------------------
 *  MAIN TEST
 *-----------------------------------------------------------------------------*/
//int main(int argc, char *argv[])
//{
//	size_t n = 10000;
//	PBar pb(n);
//
//    std::cout << "test short loop" << std::endl;
//    pb.set_description("short loop");
//    // neat to also increment the pb during the for declaration!!
//	for(size_t i=0; i <= n; ++i, ++pb) {
//		usleep(200);
//	}
//
//    pb.reset();
//    std::cout << "test short loop unknown length" << std::endl;
//    pb.set_description("short loop");
//    // neat to also increment the pb during the for declaration!!
//	for(size_t i=0; i <= n; ++i, ++pb) {
//		usleep(200);
//	}
//
//    std::cout << "\ntest long loop (but only the start)" << std::endl;
//    
//	n = 99999;
//	pb.reset(n);
//	for(size_t i=0; i <= 15; ++i, ++pb) {
//		sleep(1);
//	}
//
//    return 0;
//}
//
// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
