//
//  shutup.cpp
//  Extragas
//
//  Created by Anqi Li on 27/08/2021.
//

#include "shutup.hpp"
#include <stdio.h>
#include <unistd.h>

static char buf[20];
static int saved_stdout;

void stdout_off() {
    saved_stdout = dup(1);
    freopen("/dev/null", "w", stdout);
}

void stdout_on() {
    sprintf(buf, "/dev/fd/%d", saved_stdout);
    freopen(buf, "w", stdout);
}
