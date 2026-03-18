#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>
#include <mutex>

enum class LogLevel {
    DEBUG = 0,
    INFO = 1,
    WARN = 2,
    ERROR = 3
};

class Logger {
public:
    static void set_level(LogLevel level) {
        std::lock_guard<std::mutex> lock(get_mutex());
        get_current_level() = level;
    }

    static void debug(const std::string& msg) { log(LogLevel::DEBUG, "DEBUG: " + msg); }
    static void info(const std::string& msg)  { log(LogLevel::INFO,  "INFO:  " + msg); }
    static void warn(const std::string& msg)  { log(LogLevel::WARN,  "WARN:  " + msg); }
    static void error(const std::string& msg) { log(LogLevel::ERROR, "ERROR: " + msg); }

private:
    static LogLevel& get_current_level() {
        static LogLevel level = LogLevel::INFO;
        return level;
    }

    static std::mutex& get_mutex() {
        static std::mutex mtx;
        return mtx;
    }

    static void log(LogLevel level, const std::string& msg) {
        std::lock_guard<std::mutex> lock(get_mutex());
        if (level >= get_current_level()) {
            if (level >= LogLevel::WARN) {
                std::cerr << msg << std::endl;
            } else {
                std::cout << msg << std::endl;
            }
        }
    }
};

#endif // LOGGER_H
