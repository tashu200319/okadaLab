#!/usr/bin/env python3
import json
import time
import os

_DEBUG_LOG_PATH = "/Users/tashiroshuya/Desktop/okadaLab/.cursor/debug.log"
def _dlog(hid, loc, msg, data):
    try:
        os.makedirs(os.path.dirname(_DEBUG_LOG_PATH), exist_ok=True)
        with open(_DEBUG_LOG_PATH, "a") as f:
            f.write(json.dumps({"sessionId": "debug-session", "runId": "run1", "hypothesisId": hid, "location": loc, "message": msg, "data": data, "timestamp": round(time.time() * 1000)}) + "\n")
    except Exception as e:
        print(f"Log error: {e}")

_dlog("TEST", "test_debug:main", "Test log entry", {"test": True})
print("Test log written")
