#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define VERSION_MAJOR 1
#define VERSION_MINOR 20
#define VERSION_PATCH 0


// Helper macros to stringify values
#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x)

#define SPINWALK_VERSION STRINGIFY(VERSION_MAJOR) "." STRINGIFY(VERSION_MINOR) "." STRINGIFY(VERSION_PATCH)

#define ERR_MSG  "\033[1;31mError:\033[0m "
#define WARN_MSG "\033[1;33mWarning:\033[0m "

// #define MINI_CASE_SENSITIVE

#define GAMMA  267515315. // rad/s.T

#define SPINWALK_ASCII_ART \
" ____            _          __        __          _   _        \n" \
"/ ___|   _ __   (_)  _ __   \\ \\      / /   __ _  | | | | __    \n" \
"\\___ \\  | '_ \\  | | | '_ \\   \\ \\ /\\ / /   / _` | | | | |/ /    \n" \
" ___) | | |_) | | | | | | |   \\ V  V /   | (_| | | | |   <     \n" \
"|____/  | .__/  |_| |_| |_|    \\_/\\_/     \\__,_| |_| |_|\\_\\    \n" \
"        |_|                                                    \n\n"



#endif // DEFINITIONS_H
