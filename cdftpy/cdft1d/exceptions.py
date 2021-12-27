class ConvergenceError(Exception):
    def __init__(self, message=None):
        if message is None:
            message = "Convergence error"
        super().__init__(message)

    # def __str__(self):
    # #     print('calling str')
    #     if self.message:
    #         return 'Convergence error, {0} '.format(self.message)
    #     else:
    #         return 'Convergence error'
