# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Logging

struct TeeLogger <: AbstractLogger
    loggers::NTuple{2,AbstractLogger}
end

Logging.min_enabled_level(t::TeeLogger) = minimum(Logging.min_enabled_level.(t.loggers))
Logging.shouldlog(t::TeeLogger, level, _module, group, id) =
    any(l -> Logging.shouldlog(l, level, _module, group, id), t.loggers)
Logging.catch_exceptions(t::TeeLogger) = any(Logging.catch_exceptions, t.loggers)
function Logging.handle_message(t::TeeLogger, level, message, _module, group, id, file, line; kwargs...)
    for l in t.loggers
        if Logging.shouldlog(l, level, _module, group, id)
            Logging.handle_message(l, level, message, _module, group, id, file, line; kwargs...)
        end
    end
end

"""
    setup_logger(io_cfg; default_log_path="run.log", default_console_log=true, default_console_level="info")

Return a NamedTuple with `logger`, `logio`, and `log_path`.
"""
function setup_logger(io_cfg::AbstractDict;
    default_log_path::AbstractString="run.log",
    default_console_log::Bool=true,
    default_console_level::AbstractString="info")

    log_path = get(io_cfg, "log_path", default_log_path)
    console_log = get(io_cfg, "console_log", default_console_log)
    level_str = lowercase(String(get(io_cfg, "console_level", default_console_level)))
    console_level = level_str == "debug" ? Logging.Debug :
        level_str == "warn" ? Logging.Warn :
        level_str == "error" ? Logging.Error :
        Logging.Info

    logio = open(log_path, "a")
    file_logger = SimpleLogger(logio, Logging.Info)
    logger = console_log ?
        TeeLogger((ConsoleLogger(stderr, console_level), file_logger)) :
        file_logger

    return (; logger, logio, log_path)
end
