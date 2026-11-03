"""Deterministic playback clock; the physical solution is evaluated separately."""
from dataclasses import dataclass


@dataclass
class Playback:
    elapsed_s: float = 0.0
    running: bool = False
    last_tick: float = 0.0

    def advance(self, now, duration, rate=1.0):
        if self.running:
            self.elapsed_s = min(duration, self.elapsed_s + max(0.0, now-self.last_tick)*rate)
            if self.elapsed_s >= duration:
                self.running = False
        self.last_tick = now

    def start(self, now, duration):
        if self.elapsed_s >= duration:
            self.elapsed_s = 0.0
        self.running = True
        self.last_tick = now

    def pause(self):
        self.running = False

    def stop(self, now):
        self.elapsed_s = 0.0
        self.running = False
        self.last_tick = now
