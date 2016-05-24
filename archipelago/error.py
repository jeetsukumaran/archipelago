

class ArchipelagoException(Exception):
    pass

class FailedSimulationException(ArchipelagoException):
    pass

class InsufficientFocalAreaLineagesSimulationException(FailedSimulationException):
    pass

class TotalExtinctionException(FailedSimulationException):
    pass

class FailedToMeetFocalAreaLineageTargetException(FailedSimulationException):
    pass
