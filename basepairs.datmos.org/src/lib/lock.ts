export class AsyncLock {

    #queue: (() => void)[] = []
    #locked = false
    // #lock: Promise<void> | null = null
    // #nWaiting = 0
    #cancelOnWait: AbortController | null = null

    acquire(abort?: AbortSignal): Promise<() => void> | (() => void) {
        const runNext = () => {
            if (this.#queue.length > 0) {
                const next = this.#queue.shift()!
                next()
            }
        }
        if (!this.#locked) {
            this.#locked = true
            return () => {
                this.#locked = false
                runNext()
            }
        }

        return new Promise<() => void>((resolve, reject) => {
            const handle = () => {
                if (abort?.aborted) {
                    reject(new DOMException("The operation was aborted.", "AbortError"))
                    return;
                }
                if (this.#locked) {
                    throw new Error("Well fuck")
                }
                this.#locked = true
                resolve(() => {
                    this.#locked = false
                    runNext()
                })
            }
            this.#queue.push(handle)

            if (abort) {
                abort.addEventListener("abort", ev => {
                    const ix = this.#queue.indexOf(handle)
                    if (ix >= 0) {
                        this.#queue.splice(ix, 1)
                    }
                    reject(new DOMException("The operation was aborted.", "AbortError"))
                })
            }
        })
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
        const ac = this.#cancelOnWait ??= new AbortController()
        const release = await this.acquire(ac.signal);
        try {
            ac.signal.throwIfAborted()
            return await fn(ac.signal);
        } finally {
            release();
        }
    }

    abortRunning() {
        if (this.#cancelOnWait) this.#cancelOnWait.abort()
        this.#cancelOnWait = null
    }
}
