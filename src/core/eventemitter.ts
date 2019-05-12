type Listener<T extends Array<any>> = (...args: T) => void

export class Emitter<T extends Array<any>> {
  subscriptions: {
    [key: string]: Array<Listener<T>>
  } = {}

  emit(key: string, ...args: T) {
    const listeners = this.subscriptions[key]

    if (listeners) {
      for (var i = 0; i < listeners.length; i++) {
        listeners[i](...args)
      }
    }
  }

  on(key: string, listener: Listener<T>) {
    if (this.subscriptions[key]) {
      this.subscriptions[key].push(listener)
    } else {
      this.subscriptions[key] = [listener]
    }
  }

  off(key: string, listener: Listener<T>) {
    this.subscriptions[key] = this.subscriptions[key].filter((l) => l === listener)
  }
}
