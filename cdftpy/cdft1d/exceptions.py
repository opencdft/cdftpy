class ConvergenceError(Exception):
    def __init__(self, message=None):
        if message is None:
            message = "Convergence error"
        else:
            message = F"Convergence error: {message}"
        super().__init__(message)
