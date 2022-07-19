from abc import ABC


class NetworkStrategy(ABC):
    ...


class PriorityQueueStrategy(NetworkStrategy):
    def __init__(self) -> None:
        ...
