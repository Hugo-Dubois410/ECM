'''
Created on 5 mars 2017

@author: garzol
'''

from multiprocessing import Process, Pipe
from multiprocessing.connection import wait

def foo(w,f):
    w.send(f())
    w.close()

def parallelize(f,nProc,n):
    readers = []
    processes = []
    ret = 0
    for _ in range(nProc):
        r, w = Pipe(duplex=False)
        p = Process(target=foo, args=(w,f,))
        p.start()
        readers.append(r)
        processes.append(p)
        # We close the writable end of the pipe now to be sure that
        # p is the only process which owns a handle for it.  This
        # ensures that when p closes its handle for the writable end,
        # wait() will promptly report the readable end as being ready.
        w.close()

    while readers:
        for r in wait(readers):
            try:
                msg = r.recv()
            except EOFError:
                readers.remove(r)
            else:
                r.close()
                readers.remove(r)
                if msg>1 and msg!=n: 
                    ret = msg
                    for p in processes: p.terminate()
                    for r in readers: r.close()
                    readers = []
                    break
    return ret
