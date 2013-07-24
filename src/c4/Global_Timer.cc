//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Global_Timer.cc
 * \author Kent G. Budge
 * \date   Mon Mar 25 17:35:07 2002
 * \brief  Define methods of class Global_Timer, a POSIX standard timer.
 * \note   Copyright (C) 2002-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Timer.hh 7075 2013-04-01 22:48:15Z kellyt $
//---------------------------------------------------------------------------//

#include <iostream>

#include "Global_Timer.hh"

namespace rtt_c4
{
using namespace std;

bool Global_Timer::global_active_ = false;
map<string, Global_Timer::timer_entry> Global_Timer::active_list_;
     
//---------------------------------------------------------------------------------------//
Global_Timer::Global_Timer(char const *name)
    :
    name_(name)
{
    Require(name != NULL);

    timer_entry &entry = active_list_[name];
    active_ = entry.is_active;
    Check(entry.timer == NULL);
    entry.timer = this;

    Ensure(name == this->name());
}

//---------------------------------------------------------------------------------------//
/* static */
void Global_Timer::set_global_activity(bool const active)
{
    if (rtt_c4::node()==0)
    {
        global_active_ = active;
        
        cout << "***** Global timers are now ";
        if (active)
            cout << "ACTIVE";
        else
            cout << "INACTIVE";
        
        cout << endl;
    }
}

//---------------------------------------------------------------------------------------//
Global_Timer::~Global_Timer()
{
    if (active_ || global_active_)
    {
        cout << endl;
        cout << "Timing report for timer " << name_ << ':' << endl;
        print(cout);
        cout << endl;
    }
}

//---------------------------------------------------------------------------------------//
/*static*/
void Global_Timer::set_global_activity(set<string> const &timer_list)
{
    if (rtt_c4::node()==0)
    {
        cout << "***** Global timers selectively activated:";
        for (set<string>::const_iterator i=timer_list.begin(); i!=timer_list.end(); ++i)
        {
            string const &name = (*i);
            cout << " \"" << name << '\"';
            timer_entry &entry = active_list_[name];
            entry.is_active = true;
            if (entry.timer != NULL)
            {
                entry.timer->set_activity(true);
            }
            else
            {
                cout << " (DEFERRED)";
            }
        }
        cout << endl;
    }
}


} // end namespace rtt_c4

//---------------------------------------------------------------------------//
//                              end of c4/Global_Timer.cc
//---------------------------------------------------------------------------//