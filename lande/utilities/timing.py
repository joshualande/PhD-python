
from datetime import datetime
from dateutil import relativedelta
class StopWatch(object):
    """ Simple object for timing code"""
    def __init__(self):
        self.start_time = self.__lap = datetime.now()
  
    @property
    def lap(self):
        current_time = datetime.now()
        lap = relativedelta.relativedelta(current_time,self.__lap)
        self.__lap = current_time
        return lap

    @property
    def total(self):
        return datetime.now() - self.start_time
    
    def __str__(self):
        ret=[]
        lap=self.lap
        if lap.years > 0: ret.append('%g years' % lap.years)
        if lap.months > 0: ret.append('%g months' % lap.months)
        if lap.days > 0: ret.append('%g days' % lap.days)
        if lap.hours > 0: ret.append('%g hours' % lap.hours)
        if lap.minutes > 0: ret.append('%g minutes' % lap.minutes)
        if lap.seconds > 0: ret.append('%g seconds' % lap.seconds)

        # only add microseconds if necessary
        if len(ret) < 1: ret.append('%g microseconds' % lap.microseconds)

        return 'Lap time is '+', '.join(ret[0:2]) # only print out leading two 



