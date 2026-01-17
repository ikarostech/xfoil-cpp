#pragma once

#include <ostream>
#include <string>

class Logger {
 public:
  virtual ~Logger() = default;
  virtual void write(const std::string& message) = 0;
};

class StreamLogger : public Logger {
 public:
  explicit StreamLogger(std::ostream& out) : out_(out) {}
  void write(const std::string& message) override { out_ << message; }

 private:
  std::ostream& out_;
};
