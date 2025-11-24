import os
import sys
import logging
import datetime
from pathlib import Path
from typing import Optional, Dict, Any

class TransRateLogger:
    def __init__(self, log_dir: Path, assembly_name: str, quiet: bool = False):
        self.log_dir       = Path(log_dir)
        self.assembly_name = assembly_name
        self.quiet         = quiet
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.main_logger   = self._setup_main_logger()
        self.tool_loggers: Dict[str, logging.Logger] = {}
        
    def _setup_main_logger(self) -> logging.Logger:
        logger = logging.getLogger('transrate2_main')
        logger.setLevel(logging.INFO)
        
        logger.handlers.clear()
        
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file  = self.log_dir / f"transrate2_{self.assembly_name}_{timestamp}.log"
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        
        logger.addHandler(file_handler)

        logger.info(f"=== TransRate2 Analysis Started for {self.assembly_name} ===")
        return logger
    
    def _setup_tool_logger(self, tool_name: str) -> logging.Logger:
        logger_name = f'transrate2_{tool_name}'
        logger      = logging.getLogger(logger_name)
        logger.setLevel(logging.DEBUG)
        
        logger.handlers.clear()
        
        stdout_file = self.log_dir / f"{tool_name}_{self.assembly_name}_stdout.log"
        stderr_file = self.log_dir / f"{tool_name}_{self.assembly_name}_stderr.log"
        
        stdout_handler   = logging.FileHandler(stdout_file)
        stdout_handler.setLevel(logging.INFO)
        stdout_formatter = logging.Formatter('%(asctime)s - STDOUT - %(message)s')
        stdout_handler.setFormatter(stdout_formatter)
        
        stderr_handler   = logging.FileHandler(stderr_file)
        stderr_handler.setLevel(logging.WARNING)
        stderr_formatter = logging.Formatter('%(asctime)s - STDERR - %(message)s')
        stderr_handler.setFormatter(stderr_formatter)
        
        logger.addHandler(stdout_handler)
        logger.addHandler(stderr_handler)
        
        return logger
    
    def get_tool_logger(self, tool_name: str) -> logging.Logger:
        if tool_name not in self.tool_loggers:
            self.tool_loggers[tool_name] = self._setup_tool_logger(tool_name)
        return self.tool_loggers[tool_name]
    
    def log_progress(self, message: str, level: str = 'info'):
        if level.lower() == 'info':
            self.main_logger.info(message)
        elif level.lower() == 'warning':
            self.main_logger.warning(message)
        elif level.lower() == 'error':
            self.main_logger.error(message)
        elif level.lower() == 'debug':
            self.main_logger.debug(message)
    
    def log_tool_start(self, tool_name: str, command: list, mode: str = ""):
        logger  = self.get_tool_logger(tool_name)
        cmd_str = ' '.join(str(x) for x in command)
        
        logger.info(f"=== {tool_name.upper()} {mode} STARTED ===")
        logger.info(f"Command: {cmd_str}")
        logger.info(f"Working directory: {os.getcwd()}")
        
        self.log_progress(f"Starting {tool_name} {mode}")
    
    def log_tool_output(self, tool_name: str, stdout: str, stderr: str, returncode: int):
        logger = self.get_tool_logger(tool_name)
        
        if stdout:
            for line in stdout.strip().split('\n'):
                if line:
                    logger.info(line)
        
        if stderr:
            for line in stderr.strip().split('\n'):
                if line:
                    logger.warning(line)
        
        logger.info(f"=== Process completed with return code: {returncode} ===")
        
        if returncode == 0:
            self.log_progress(f"{tool_name} completed successfully")
        else:
            self.log_progress(f"{tool_name} failed with return code {returncode}", 'error')
    
    def log_stage_start(self, stage_name: str, details: Optional[Dict[str, Any]] = None):
        message = f"=== Starting {stage_name} ==="
        if details:
            detail_str  = ", ".join(f"{k}: {v}" for k, v in details.items())
            message    += f" ({detail_str})"
        self.log_progress(message)
    
    def log_stage_complete(self, stage_name: str, metrics: Optional[Dict[str, Any]] = None):
        message = f"=== Completed {stage_name} ==="
        if metrics:
            metric_str  = ", ".join(f"{k}: {v}" for k, v in metrics.items())
            message    += f" - Metrics: {metric_str}"
        self.log_progress(message)
    
    def log_file_operation(self, operation: str, file_path: str, success: bool = True):
        status = "SUCCESS" if success else "FAILED"
        self.log_progress(f"File {operation}: {file_path} - {status}")
    
    def finalize(self):
        self.log_progress("=== TransRate2 Analysis Completed ===")
        
        for handler in self.main_logger.handlers[:]:
            handler.close()
            self.main_logger.removeHandler(handler)
        
        for tool_logger in self.tool_loggers.values():
            for handler in tool_logger.handlers[:]:
                handler.close()
                tool_logger.removeHandler(handler)

class LoggingSubprocess:
    def __init__(self, logger: TransRateLogger, tool_name: str):
        self.logger    = logger
        self.tool_name = tool_name
    
    def run_with_logging(self, command: list, mode: str = "") -> tuple:
        import subprocess
        
        self.logger.log_tool_start(self.tool_name, command, mode)
        
        try:
            process = subprocess.Popen(
                command,
                stdout     = subprocess.PIPE,
                stderr     = subprocess.PIPE,
                preexec_fn = os.setsid,
                shell      = False,
                text       = True
            )
            
            stdout, stderr = process.communicate()
            returncode     = process.returncode
            
            self.logger.log_tool_output(self.tool_name, stdout, stderr, returncode)
            
            return returncode, stdout, stderr
            
        except Exception as e:
            error_msg = f"Failed to execute {self.tool_name}: {str(e)}"
            self.logger.log_progress(error_msg, 'error')
            return 1, "", error_msg
