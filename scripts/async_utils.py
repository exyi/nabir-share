
class MockPool:
    class ApplyResult:
        def __init__(self, x) -> None:
            self.x = x
        def get(self):
            return self.x
    def apply_async(self, f, args, callback=None):
        r = f(*args)
        if callback is not None:
            callback(r)
        return MockPool.ApplyResult(r)
