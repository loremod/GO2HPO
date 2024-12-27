class LogManager:
    def __init__(self, is_active=True):
        self.is_silent = not is_active
    def activate(self):
        self.is_silent = False
    def deactivate(self):
        self.is_silent = True
    def toggle(self):
        self.is_silent = not self.is_silent
    def log(self, *args):
        if not self.is_silent:
            print(*args)