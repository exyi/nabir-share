export class AsyncLock {

    #lock: Promise<void> | null = null
    #nWaiting = 0
    #cancelOnWait: AbortController | null = null

    async acquire(): Promise<() => void> {
        this.#nWaiting += 1;
        while (this.#lock !== null) {
            this.abortRunning()
            await this.#lock;
        }

        if (this.#lock !== null) throw 'wtf'

        let resolve: null | (() => void) = null
        this.#lock = new Promise(r => resolve = r);
        
        this.#nWaiting -= 1;

        if (!resolve) throw 'wtf'
        return () => {
            this.#lock = null;
            resolve();
        }
    }

    async withLock<T>(fn: () => Promise<T>): Promise<T> {
        const release = await this.acquire();
        try {
          return await fn();
        } finally {
          release();
        }
    }

    async withCancellableLock<T>(fn: (signal: AbortSignal) => Promise<T>): Promise<T> {
        const release = await this.acquire();
        try {
            const controller = new AbortController()
            this.#cancelOnWait = controller
            return await fn(controller.signal);
        } finally {
            this.#cancelOnWait = null
            release();
        }
    }

    abortRunning() {
        if (this.#cancelOnWait) this.#cancelOnWait.abort()
        return this.#lock
    }
}
