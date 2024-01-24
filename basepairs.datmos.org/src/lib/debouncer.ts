export class AsyncDebouncer {
    #timeout: number
    #animationFrame: boolean
    constructor(timeout: number, animationFrame = true) {
        this.#timeout = timeout
        this.#animationFrame = animationFrame
    }

    #timerId: any = null
    #animationFrameId: any = null

    #resolve: (() => void) | null = null
    #promise: Promise<void> | null = null;
    #resultPromise: Promise<void> | null = null;


    resetTimer() {
        if (!this.#promise) {
            this.#promise = new Promise(resolve => {
                this.#resolve = resolve
            })
        }
        const resolve = () => {
            if (this.#resolve) {
                this.#resolve()
                this.#resolve = null
            }
            this.#promise = null
            this.#resultPromise = null
        }

        if (this.#animationFrameId) {
            cancelAnimationFrame(this.#animationFrameId)
            this.#animationFrameId = null
        }

        if (this.#timerId) {
            clearTimeout(this.#timerId)
        }

        this.#timerId = setTimeout(() => {
            this.#timerId = null
            if (this.#animationFrame) {
                this.#animationFrameId = requestAnimationFrame(() => {
                    this.#animationFrameId = null
                    resolve()
                })
            } else {
                resolve()
            }
        }, this.#timeout)
    }

    debounceCheckpoint(): Promise<void> {
        this.resetTimer()
        return this.#promise
    }

    debounce(fn: () => void): Promise<void> {
        this.resetTimer()

        if (this.#resultPromise) {
            return this.#resultPromise
        }

        return this.#resultPromise = this.#promise.then(fn)
    }
}
