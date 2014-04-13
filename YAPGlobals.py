"""Variables to be used globally across YAP modules"""

## leave batch job scripts and standard output (this can be controled
## through command line as well
debug_grid_tasks = False

## this will cause exception when StepDefault.find() does not find
## anything, but it turns out to be a routine situation, so do not
## use this yet
find_require = False

## use this if you need to run interactive debugger (e.g. pdb) on
## steps. Steps will inherit from dummy_threading.Thread, and thus
## run in the main thread sequentially. Note that batch job monitoring
## and logging will be still real threads. Currently, both will run
## forever with this switch. Just kill the python process once you
## are done debugging.
step_dummy_thread = False

## maximum number of concurrent Step threads (through semaphore count)
step_threads_max = 100 #500

## maximum number of concurrently opened files (through semaphore count)
open_files_max = 200 #400

