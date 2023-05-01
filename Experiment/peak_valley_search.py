# In[pkgs]:
    
from collections import Counter
from scipy.signal import detrend
import numpy as np
from enum import IntEnum

# In[peak valley search]

class Event(IntEnum):
    PEAK = 0
    VALLEY = 1


class EventState:
    def __init__(self, event, close_prices, best_price, best_idx=0):
        self.event = event
        self._close_prices = close_prices
        self.best_price = best_price
        self.best_idx = best_idx
        self._event_indices = {Event.PEAK: [], Event.VALLEY: []}

    @property
    def event_inices(self):
        return tuple(self._event_indices.values())

    def switch(self):
        """Switch event state (name)"""
        self.event ^= 1

    def set_best_rate(self, best_idx):
        self.best_price = self._close_prices[best_idx]
        self.best_idx = best_idx

    def add_best_index(self):
        self._event_indices[self.event].append(self.best_idx)


def get_peak_valley(arr, threshold, window_size, overlap, req_angles):
    # validate params
    window_size = round(window_size)
    req_angles = max(round(req_angles), 1)
    window_step = max(round(window_size * (1 - overlap)), 1)

    # get all points that classify as a peak/valley
    peak_counts, valley_counts = Counter(), Counter()
    arr_size = len(arr)

    for i in range(0, arr_size, window_step):
        flattened = detrend(arr[i:i + window_size])
        std, avg = np.std(flattened), np.mean(flattened)
        lower_b = avg - std * threshold
        upper_b = avg + std * threshold

        for idx, val in enumerate(flattened):
            if val < lower_b:
                valley_counts[idx + i] += 1
            elif val > upper_b:
                peak_counts[idx + i] += 1

    # discard points that have counts below the threshold
    pk_inds = [i for i, c in peak_counts.items() if c >= req_angles]
    vly_inds = [i for i, c in valley_counts.items() if c >= req_angles]

    # initialize iterator to find to best peak/valley for consecutive detections
    if len(pk_inds) == 0 or len(vly_inds) == 0:
        return pk_inds, vly_inds

    if pk_inds[0] < vly_inds[0]:
        curr_event, best_price = Event.PEAK, close_prices[pk_inds[0]]
    else:
        curr_event, best_price = Event.VALLEY, close_prices[vly_inds[0]]

    event_state = EventState(curr_event, close_prices=close_prices, best_price=best_price)
    event_inds = sorted(pk_inds + vly_inds)
    peak_ids_set = set(pk_inds)

    # iterate through points and only carry forward the index
    # that has the highest or lowest value from the current group
    for x in event_inds:
        in_peak = x in peak_ids_set
        is_peak = event_state.event == Event.PEAK
        is_valley = event_state.event == Event.VALLEY

        if (in_peak and is_valley) or (not in_peak and is_peak):
            event_state.add_best_index()
            event_state.switch()
            event_state.set_best_rate(best_idx=x)
            continue

        if (in_peak and is_peak and close_prices[x] > event_state.best_price) or \
                (not in_peak and is_valley and close_prices[x] < event_state.best_price):
            event_state.set_best_rate(x)

    event_state.add_best_index()

    return event_state.event_inices