//
// Created by Mert Biyikli on 25.06.26.
//

#ifndef OBLIVIOUSROUTING_ERRORS_H
#define OBLIVIOUSROUTING_ERRORS_H

#include <expected>
#include <string>

enum class ErrorCode {
    FileNotFound,
    FormatNotFound,
    InvalidGraph,
    InputNotFound,
    InvalidDemand,
    InvalidSolver,
    InvalidRouting,
    SolverFailed,
    RuntimeError,
    LogicError,
    InvalidArgument,
    UnknownException,
    NumericalFailure,
    UnsupportedSolver
};

struct Error {
    ErrorCode code;
    std::string message;
};

template <typename T>
using Result = std::expected<T, Error>;

inline auto makeErrorMessage(const ErrorCode& error, const std::string& message) {
    return std::unexpected(Error{error, message});
}

// Convert a std::exception to the project's Error type with a chosen ErrorCode.
inline auto fromStdException(const std::exception& e, ErrorCode code = ErrorCode::UnknownException) {
    return makeErrorMessage(code, std::string(e.what()));
}

// Fallback for non-std exceptions caught via catch(...)
inline auto fromUnknownException(ErrorCode code = ErrorCode::UnknownException) {
    return makeErrorMessage(code, "Unknown exception");
}

template <typename T>
inline auto getError(Result<T>& error) {
    return std::unexpected(error.error());
}

#endif //OBLIVIOUSROUTING_ERRORS_H