#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>

class Logger {
 public:
  virtual ~Logger() = default;
  virtual void write(const std::string& message) = 0;
  static Logger& instance();
};

class NullLogger : public Logger {
 public:
  void write(const std::string& message) override {
    (void)message;
  }
};

class StreamLogger : public Logger {
 public:
  explicit StreamLogger(std::ostream& out) : out_(out) {}
  void write(const std::string& message) override { out_ << message; }

 private:
  std::ostream& out_;
};

class FileLogger : public Logger {
 public:
  explicit FileLogger(const std::string& path)
      : out_(path, std::ios::app) {}
  bool isOpen() const { return out_.is_open(); }
  void write(const std::string& message) override { out_ << message; }

 private:
  std::ofstream out_;
};

inline Logger& Logger::instance() {
  static std::unique_ptr<Logger> logger = []() -> std::unique_ptr<Logger> {
    const char* env_value = std::getenv("XFOIL_LOG");
    if (!env_value || env_value[0] == '\0') {
      return std::make_unique<NullLogger>();
    }
    const std::string destination(env_value);
    if (destination == "stdout") {
      return std::make_unique<StreamLogger>(std::cout);
    }
    if (destination == "stderr") {
      return std::make_unique<StreamLogger>(std::cerr);
    }
    auto file_logger = std::make_unique<FileLogger>(destination);
    if (!file_logger->isOpen()) {
      return std::make_unique<StreamLogger>(std::cerr);
    }
    return file_logger;
  }();
  return *logger;
}
