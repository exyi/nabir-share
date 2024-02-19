
class MockPool:
    class ApplyResult:
        def __init__(self, x) -> None:
            self.x = x
        def get(self):
            return self.x
    def apply_async(self, f, args):
        return MockPool.ApplyResult(f(*args))
