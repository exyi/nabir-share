export class Lazy<T> {
    #value: T | null = null
    #fn: () => T | null
    constructor(fn: () => T) {
        if (fn == null) throw new Error("fn must not be null")
        this.#fn = fn
    }
    get isComputed() {
        return this.#fn == null
    }
    get value() {
        if (this.#fn == null)
            return this.#value
        else {
            const fn = this.#fn
            this.#fn = null
            this.#value = fn()
            return this.#value
        }
    }
}
