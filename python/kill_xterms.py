import os
import sys
import signal

# process utilities
import psutil as PS


kill_them = True

# find running processes

p_name = 'xterm'

scr_name = 'run_all_d'

ip = []

procs = [p for p in PS.process_iter()]
for i,p in enumerate(procs):
    if (p.name() == p_name) and (p.cmdline()[-1].find(scr_name) >= 0):
        if kill_them:
            print(' killing : ', p.cmdline())
            p.kill()
        else:
            print(' found : ', p.cmdline())





