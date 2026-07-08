#pragma once
#include <string>
#include <unordered_map>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include <Rcpp.h>

namespace nm_pred {

struct ExprEnv {
  std::unordered_map<std::string, double> vars;
  Rcpp::NumericVector theta;
  Rcpp::NumericVector eta;
};

class ExprParser {
 public:
  explicit ExprParser(const std::string& s, ExprEnv* env)
    : src_(s), pos_(0), env_(env) {}

  static double eval_with_env(const std::string& expr, ExprEnv& env) {
    ExprParser p(expr, &env);
    double v = p.parse_expr();
    p.skip_ws();
    if (p.pos_ < p.src_.size()) {
      throw std::runtime_error("Unexpected trailing characters in expression.");
    }
    return v;
  }

  static bool is_parseable(const std::string& expr) {
    try {
      Rcpp::NumericVector th(10), et(10);
      for (int i = 0; i < 10; ++i) { th[i] = 1.0; et[i] = 0.0; }
      ExprEnv env;
      env.theta = th;
      env.eta = et;
      (void)eval_with_env(expr, env);
      return true;
    } catch (...) {
      return false;
    }
  }

 private:
  std::string src_;
  size_t pos_;
  ExprEnv* env_;

  void skip_ws() {
    while (pos_ < src_.size() && std::isspace(static_cast<unsigned char>(src_[pos_]))) ++pos_;
  }

  bool match(const std::string& s) {
    skip_ws();
    if (src_.compare(pos_, s.size(), s) == 0) {
      pos_ += s.size();
      return true;
    }
    return false;
  }

  double parse_expr() {
    double v = parse_term();
    while (true) {
      skip_ws();
      if (match("+")) v += parse_term();
      else if (match("-")) v -= parse_term();
      else break;
    }
    return v;
  }

  double parse_term() {
    double v = parse_unary();
    while (true) {
      skip_ws();
      if (match("*")) v *= parse_unary();
      else if (match("/")) v /= parse_unary();
      else break;
    }
    return v;
  }

  double parse_unary() {
    skip_ws();
    if (match("+")) return parse_unary();
    if (match("-")) return -parse_unary();
    return parse_power();
  }

  double parse_power() {
    double v = parse_atom();
    skip_ws();
    while (true) {
      if (match("**")) {
        v = std::pow(v, parse_unary());
      } else if (match("^")) {
        v = std::pow(v, parse_unary());
      } else {
        break;
      }
      skip_ws();
    }
    return v;
  }

  std::string parse_ident() {
    skip_ws();
    size_t start = pos_;
    if (pos_ >= src_.size() || !(std::isalpha(static_cast<unsigned char>(src_[pos_])) || src_[pos_] == '_')) {
      throw std::runtime_error("Expected identifier.");
    }
    ++pos_;
    while (pos_ < src_.size() &&
           (std::isalnum(static_cast<unsigned char>(src_[pos_])) || src_[pos_] == '_')) {
      ++pos_;
    }
    return src_.substr(start, pos_ - start);
  }

  double parse_atom() {
    skip_ws();
    if (match("(")) {
      double v = parse_expr();
      if (!match(")")) throw std::runtime_error("Expected ')'");
      return v;
    }
    if (std::isdigit(static_cast<unsigned char>(src_[pos_])) ||
        (src_[pos_] == '.' && pos_ + 1 < src_.size() && std::isdigit(static_cast<unsigned char>(src_[pos_ + 1])))) {
      size_t start = pos_;
      while (pos_ < src_.size()) {
        char c = src_[pos_];
        if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '.' ||
              c == 'e' || c == 'E' || c == '+' || c == '-')) break;
        ++pos_;
      }
      return std::stod(src_.substr(start, pos_ - start));
    }
    std::string id = parse_ident();
    std::string up = id;
    for (char& c : up) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (up == "EXP" || up == "LOG" || up == "SQRT" || up == "ABS") {
      if (!match("(")) throw std::runtime_error("Expected '(' after function.");
      double arg = parse_expr();
      if (!match(")")) throw std::runtime_error("Expected ')'");
      if (up == "EXP") return std::exp(arg);
      if (up == "LOG") return std::log(arg);
      if (up == "SQRT") return std::sqrt(arg);
      return std::fabs(arg);
    }
    if (up == "THETA") {
      if (!match("(")) throw std::runtime_error("Expected '(' after THETA.");
      double idx = parse_expr();
      if (!match(")")) throw std::runtime_error("Expected ')'");
      int i = static_cast<int>(std::round(idx));
      if (i < 1 || i > env_->theta.size()) return 0.0;
      return env_->theta[i - 1];
    }
    if (up == "ETA") {
      if (!match("(")) throw std::runtime_error("Expected '(' after ETA.");
      double idx = parse_expr();
      if (!match(")")) throw std::runtime_error("Expected ')'");
      int i = static_cast<int>(std::round(idx));
      if (i < 1 || i > env_->eta.size()) return 0.0;
      return env_->eta[i - 1];
    }
    auto it = env_->vars.find(up);
    if (it != env_->vars.end()) return it->second;
    throw std::runtime_error("Unknown identifier: " + id);
  }
};

inline std::string pred_rhs(const std::string& line) {
  size_t eq = line.find('=');
  if (eq == std::string::npos) {
    eq = line.find("<-");
    if (eq == std::string::npos) return line;
    std::string rhs = line.substr(eq + 2);
    size_t a = rhs.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = rhs.find_last_not_of(" \t\r\n");
    return rhs.substr(a, b - a + 1);
  }
  std::string rhs = line.substr(eq + 1);
  size_t a = rhs.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  size_t b = rhs.find_last_not_of(" \t\r\n");
  return rhs.substr(a, b - a + 1);
}

inline std::string pred_lhs(const std::string& line) {
  size_t eq = line.find("<-");
  if (eq != std::string::npos) {
    std::string lhs = line.substr(0, eq);
    size_t a = lhs.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = lhs.find_last_not_of(" \t\r\n");
    return lhs.substr(a, b - a + 1);
  }
  eq = line.find('=');
  if (eq == std::string::npos) return line;
  std::string lhs = line.substr(0, eq);
  size_t a = lhs.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  size_t b = lhs.find_last_not_of(" \t\r\n");
  return lhs.substr(a, b - a + 1);
}

inline std::string upper_copy(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

}  // namespace nm_pred
