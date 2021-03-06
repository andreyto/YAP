########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## A library of "steps" or program wrappers to construct pipelines
## Pipeline steps orchestration, grid management and output handling.
#################################################

import YAPGlobals
import sys, tempfile, shlex, glob, os, stat, hashlib, time, datetime, re, curses
import shutil
import threading
from threading import *
import dummy_threading
import subprocess
from subprocess import *
from MothurCommandInfoWrapper import *
from collections import defaultdict
from collections import deque
from random import *
from Queue import *

import smtplib
from email.mime.text import MIMEText
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders

import traceback
import pdb
##threading redefines enumerate() with no arguments. as a kludge, we drop it here
globals().pop('enumerate',None)

_author="Sebastian Szpakowski"
_date="2012/09/20"
_version="Version 2"


def rmrf(targ_path):
    if os.path.islink(targ_path):
        os.remove(targ_path)
    elif os.path.exists(targ_path):
        if os.path.isdir(targ_path):
            shutil.rmtree(targ_path)
        else:
            try:
                os.remove(targ_path)
            except:
                pass

def clean_flag_file(fname):
    rmrf(fname)

check_flag_str_OK = "OK"
check_flag_str_Fail = "Fail"

def check_flag_file(fname):
    ## try multiple times while the content of flag file
    ## is inconclusive or file does not exist, in order to
    ## let shared file system view get updated
    n_max_tries = 4
    sec_wait = 5
    i_try = 1
    while True:
        try:
            if os.path.isfile(fname):
                with open(fname,"r") as inp:
                    cont = inp.read()
                    cont = cont.strip() 
                    if cont == check_flag_str_OK:
                        return True
                    elif cont == check_flag_str_Fail:
                        return False
            if i_try >= n_max_tries:
                return False
            time.sleep(sec_wait)
            i_try += 1
        except:
            pass
    return False

         


def pseudo_shuffle(strings,skip=0):
    """Shuffle strings deterministically (by hash).
    Keep the original position of the first 'skip' elements"""
    it = iter(strings)
    for i in xrange(skip):
        yield next(it)
    pairs = [ (hashlib.md5("{}.{}".format(t[1],t[0])).hexdigest(),t[1]) for t in enumerate(it) ]
    for t in sorted(pairs):
        yield t[1]

#################################################
##      Classes
##


class StoppableThread(threading.Thread):
    """Thread class with a stop() method. The thread itself has to check
    regularly for the stopped() condition."""

    def __init__(self):
        super(StoppableThread, self).__init__()
        self._stop = threading.Event()

    def stop(self):
        self._stop.set()

    def stopped(self):
        return self._stop.isSet()

class StoppableDummyThread(dummy_threading.Thread):
    """Dummy Thread class with a stop() method. The thread itself has to check
    regularly for the stopped() condition."""

    def __init__(self):
        super(StoppableThread, self).__init__()
        self._stop = dummy_threading.Event()

    def stop(self):
        self._stop.set()

    def stopped(self):
        return self._stop.isSet()

class   ReportingThread(StoppableThread):
    def __init__(self):
        super(ReportingThread, self).__init__()

    def run(self):
        #attempting to catch threads dying
        #stack trace should be printed anyway
        #for non-daemon threads, so this will
        #probably give us no benefit, just trying
        try:
            self.do_run()
        except:
            traceback.print_exc()
            raise

class   ReportingDummyThread(StoppableDummyThread):
    def __init__(self):
        super(ReportingDummyThread, self).__init__()

    def run(self):
        #attempting to catch threads dying
        #stack trace should be printed anyway
        #for non-daemon threads, so this will
        #probably give us no benefit, just trying
        try:
            self.do_run()
        except:
            traceback.print_exc()
            raise


class   BufferedOutputHandler(ReportingThread):
    def __init__(self, usecurses=False):
        super(BufferedOutputHandler, self).__init__()
        self.shutdown=False
        self.cache = deque()
        self.registered=0
    
        self.ids = list()
        self.wrap = 140
        
        self.starttime = time.time()
        
        #### init log
        try:
            self.otptfile = open("logfile.txt", 'a')
            self.toPrint("-----", "GLOBAL", "Appending to a logfile.txt...")
        except:
            self.otptfile = open("logfile.txt", 'w')
            self.toPrint("-----", "GLOBAL", "Creating a new logfile.txt...")
        
        command = " ".join(sys.argv)
        self.otptfile.write("command: %s\n" % command)
        
        #### init output (curses)
        
        
        self.usecurses = usecurses
        
        if (self.usecurses):
            
            self.stdscr=curses.initscr()
            curses.savetty()
            curses.noecho()
            curses.cbreak()
            curses.curs_set(0)
            
            self.textbuffer= list()
            self.stself.stdscr.refresh()
            self.cursestrackbuffer = 100
            self.scrollpad = curses.newpad(self.cursestrackbuffer*2, self.wrap*2)
            self.spacerpad = curses.newpad(1,1000)
            self.updatepad = curses.newpad(10,1000)
            self.rows, self.cols = self.stdscr.getmaxyx()
            
        else:
            self.stdscr=None
            
        self.start()
        
    def do_run(self):
        self.toPrint("-----", "GLOBAL", "Setting up the pipeline...")
        self.flush()
        
        time.sleep(5)
        
        while YAPGlobals.step_dummy_thread or (activeCount()>3 or self.registered>0 or len(self.cache) > 0):
            self.flush()
            time.sleep(1)
            if self.stopped():
                break
        
        self.flush()
        endtime = time.time()
        text =  "+%s [fin]" % (str(datetime.timedelta(seconds=round(endtime-self.starttime,0))).rjust(17)) 
        self.toPrint("-----", "GLOBAL", text)   
        
        command = "%spython %straverser.py" % (binpath, scriptspath)
        p = Popen(shlex.split(command), stdout = PIPE, stderr = PIPE, close_fds=True)
        dot, err = p.communicate()
        
        with open("workflow.dot", "w") as x:
            x.write(dot)
            x.write("\n")
        #DEBUG:
        skipDot = True #dot was getting into endless loop
        if not skipDot:
            for format in ["svg", "svgz", "png", "pdf"]:
                command = "dot -T%s -o workflow.%s" % (format, format) 
                p = Popen(shlex.split(command), stdin = PIPE, stdout = PIPE, stderr = PIPE, close_fds=True)
                out, err = p.communicate(dot)
                
            
            self.toPrint("-----", "GLOBAL", "Check out workflow.{svg,png,jpg} for an overview of what happened.")
        else:
            self.toPrint("-----", "GLOBAL", "Skipping call to dot graphics generation")
        self.flush()
        self.otptfile.close()
        self.closeDisplay()
        self.mailLog()
        
    def register(self, id):
        
        if(id in set(self.ids)):
            msg = "CRITICAL: Attempt to register duplicate Step ID: {}".format(id)
            self.toPrint("-----", "GLOBAL", msg)
            #this is called in the main thread
            raise ValueError(msg)

        self.registered+=1
        self.ids.append(id)
    
    def deregister(self):
        self.registered-=1
    
    def collapseIDs(self, text ):
        for id in self.ids:
            if len(id)>5:
                text = re.sub(id, "[{0}~]".format(id[:5]), text)
        return (text)   
    
    def flush(self):
            
        while len(self.cache) > 0:
            id, name, line, date = self.cache.popleft()
            tag = "[{2}] [{0}] {1:<20} > ".format( id, name, time.asctime(date) ) 
            line = "{0!s}".format(line)
            #line = self.collapseIDs(line)
            
            
            otpt = "{0}{1}".format(tag, line[:self.wrap])
            self.otptfile.write("{0}{1}\n".format(tag, line))
            
            line = line[self.wrap:]
            self.outputScroll(otpt)
            
            
            while len(line)>=self.wrap:
                otpt = "{0}\t{1}".format(tag, line[:self.wrap])
                line = line[self.wrap:]
                self.outputScroll(otpt) 
                
            if len(line)>0: 
                otpt = "{0:<30}\t\t{1}".format("", line)
                line = line
                self.outputScroll(otpt) 
        
        self.redrawScreen()                     
    def mailLog(self):
        log = loadLines("logfile.txt")
        log.reverse()
        
        paths = os.getcwd()
        paths = "%s/" % (paths)
        dirs = glob.glob("*OUTPUT*")
        dirs.sort()
        
        for d in dirs:
            paths = "%s\n\t%s/*" % (paths, d)
    
        header = "Hi,\nYAP has just finished. Most, if not all, of your data should be in:\n\n%s\n\n-see the log below just to make sure...\nThe attached work-flow graph can be opened in your browser.\nYours,\n\n~YAP"  % (paths)    
        log = "".join(log)
        msgtext = "%s\n\n<LOG>\n\n%s\n</LOG>\n\n" % (header, log)
        
        try:
            me = __email__
            toaddr = [me]
            
            
            msg = MIMEMultipart()
            msg['To'] = COMMASPACE.join(toaddr)
            msg['Date'] = formatdate(localtime=True)
            msg['Subject'] = '[AUTOMATED] YAP is done.' 
            
            
            if me != __admin__:
                ccaddr = [__admin__]
                msg['BCC'] = COMMASPACE.join(ccaddr)
                toaddr = toaddr + ccaddr
            
            
            msg.attach(MIMEText(msgtext))
            
            files = ["workflow.pdf"]
            for f in files:
                try:
                    part = MIMEBase('application', "octet-stream")
                    part.set_payload( open(f,"rb").read() )
                    Encoders.encode_base64(part)
                    part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(f))
                    msg.attach(part)
                except:
                    pass
                          
            s = smtplib.SMTP('mail.jcvi.org')
            s.sendmail(me, toaddr , msg.as_string())
            s.quit()
        except:
            pass    
        
    def redrawScreen(self):
        try:
            y,x = self.stdscr.getmaxyx()
            ### enough screen to print:
            if y>20 and x>20:
                
                if len(self.textbuffer) < (y-10):
                    self.scrollpad.refresh(0, 0, 0, 0, y-10, x-5)
                else:
                    self.scrollpad.refresh(self.cursestrackbuffer-y+10 , 0, 0, 0, y-10, x-5)    
                
                self.updatepad.refresh(0, 0, y-8, 10 , y-3, x-5)
                
            ### when screen too small
            else:   
                self.scrollpad.refresh(0,0,0,0,0,0)
                self.updatepad.refresh(0,0,0,0,0,0)
        except:
            self.closeDisplay()
            self.usecurses=False        
#                                                       
    def toPrint(self, id, name, line, date=None):
        if date is None:
            date = time.localtime()
        self.cache.append((id, name, line, date))
    
    def outputScroll(self, k):  
        if self.usecurses:
            self.textbuffer.append("%s\n" %(k))
            self.scrollpad.clear()
            for k in self.textbuffer[-self.cursestrackbuffer:]:
                self.scrollpad.addstr(k)        
        else:
            print k 
        
    def outputUpdate(self,k):
        if self.usecurses:
            self.updatepad.clear()
            for k in k.strip().split("\n"):
                self.updatepad.addstr("%s\n" % k)
                
            
    def closeDisplay(self): 
        if self.usecurses:
            
            self.stdscr.clear()
            self.stdscr.refresh()
             
            curses.curs_set(1)
            curses.nocbreak()
            curses.echo()
            curses.resetty()
            curses.endwin()
            
class   TaskQueueStatus(ReportingThread):
    def __init__(self, update=1, maxnodes=10):
        if YAPGlobals.step_dummy_thread:
            self.quiet = True
        else:
            self.quiet = False

        ReportingThread.__init__(self)
        self.active=True
        
        self.maxnodes = maxnodes 
        self.available = self.maxnodes
        
        self.update = update    
            
        #### queue of grid jobs to run
        self.scheduled = Queue()
        #### to keep track of things popped off the queue
        self.processing = dict()
        
        #### inventory of what ran
        #### tuple (jid, status) indexed by command
        #### status: new/running/done/remove
        #### new         upon registering
        #### running     when submitted to the grid
        #### done        when completed
        
        self.registered = dict()
                
        #### inventory of completed jobs        
        self.bestqueue = "default.q"
        self.pollqueues()
        
        self.running=0
        self.stats=dict()
        
        self.previous =""

        ## All task submitters with wait for this condition,
        ## that will be signalled in task.setCompleted method.
        ## All waited tasks will get notified and check for their isCompleted status (that will be
        ## serialized because lock has to be acquired by Condition.wait()).
        ## The task that was set as completed will exit the wait loop.
        ## A much more straightforward use of Event associated with every
        ## Task would however consumed one handle per event object, probably
        ## leading to thread resource errors that we have seen before.
        self.any_task_completed = threading.Condition()
        
        self.start()
            
    def do_run(self):
        BOH.toPrint("-----","BATCH","Setting up the grid...")
        time.sleep(5)
        while not self.stopped():
        #while YAPGlobals.step_dummy_thread or (activeCount()>3 or self.running>0 or self.scheduled.qsize()>0):

            self.pollfinished()
            self.pollqueues()
            self.pollrunning()
            self.dispatch()
            self.cleanup()
            if not self.quiet:
                BOH.toPrint("-----","BATCH","{}".format(self))
            
            time.sleep(self.update)
            
        BOH.toPrint("-----","BATCH","{}\nGrid Offline.".format(self))
        
        print self  
        print "Queue status shutting down."
        
    def cleanup(self):
        toremove = set()
        for key, tup in self.registered.items():
            id, status = tup
            if status == "remove":
                toremove.add(key)
        for key in toremove:
            del self.registered[key]
    
    def flagRemoval(self, task):
        id, status = self.registered[task.getUniqueID()]
        if status =="done":
             self.registered[task.getUniqueID()] = [id, "remove"]
        else:
            print "cannot flag yet:", id, status
            
    def pollfinished(self):
                    
#       donejobs = set() 
#       
#       ### only 100 recent jobs shown, which could be a problem ;-)   
#       p = Popen(shlex.split("qstat -s z"), stdout=PIPE, stderr=PIPE, close_fds=True)
#       
#       out,err = p.communicate()
#
#       lines = out.split("\n")
#       tmp = set()
#       if len(lines)>2:
#           for line in lines[2:]:
#               line = line.strip().split()
#               if len(line)>0:
#                   donejobs.add(line[0])       
#       
        #if len(donejobs)>0:
        for key, tup in self.registered.items():
            id, status = tup
            #if (status == "running") and (id in donejobs):
            if (status == "running") and (self.isJobDone(id)):
                tmp = self.registered[key][1]= "done"
                self.processing[key].setCompleted()
                self.available += 1
                del self.processing[key]
                
    def isJobDone(self, jid):
        if jid == -1:
            BOH.toPrint("-----","BATCH","Impossible job id for qstat {}, marking as done...".format(jid))
            return True
        err = ""
        for i_try in range(3):
            time.sleep(2**(i_try+1)-2)
            p = Popen(shlex.split("qstat -j %s" % jid), stdout=PIPE, stderr=PIPE, close_fds=True)
            out,err = p.communicate()       
            if err.find("jobs do not exist")>-1:
                return True
            elif p.returncode == 0:
                break
            if not self.quiet:
                BOH.toPrint("-----","BATCH","qstat error {}, trying again...".format(err))
        else:
            if not self.quiet:
                BOH.toPrint("-----","BATCH","isJobDone() multiple qstat errors {}, giving up.".format(err))
        return False
           
    def pollqueues(self):

        qstat_ok = False
        command="qstat -g c" 
        err = ""
        for i_try in range(3):
            time.sleep(2**(i_try+1)-2)
            p = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, close_fds=True )
            out,err = p.communicate()
            #no point to continue if host configured to never submit
            assert err.find("neither submit nor admin host")==-1
            if p.returncode == 0:
                qstat_ok = True
                break
            if not self.quiet:
                BOH.toPrint("-----","BATCH","qstat error {}, trying again...".format(err))
        else:
            if not self.quiet:
                BOH.toPrint("-----","BATCH","pollqueues() multiple qstat errors {}, giving up.".format(err))
        
        if qstat_ok:
            #else we just keep the previous attribute values
            
            queues = defaultdict(float)
            out = out.strip().split("\n")
            
            fullqueues = set()
            
            #cache queue information
            for q in out[2:]:
                queue, cqload, used, res, avail, total, acds, cdsu = q.split()
                avail = float(avail)
                total = float(total)
                if total>0:
                    queues[queue] = avail
                
                if avail==0:
                    fullqueues.add(queue)   
            
            
            # determine which queue is the best
            
            #for k in ("default.q", "medium.q", "fast.q", "himem.q"):
            #for k in ("fast.q", "medium.q", "default.q"):
            #for k in ("himem.q", "medium.q", "default.q"):
                                
            if ("medium.q" in fullqueues) and ("default.q" in fullqueues) and "himem" in queues.keys() :
                if queues["himem.q"]>0:
                    self.bestqueue = "himem.q"
                else:
                    self.bestqueue = "medium.q"
            else:
                for k in ("medium.q", "default.q"):
                    if queues[k] >= queues[self.bestqueue]: 
                        self.bestqueue = k          
            if YAPGlobals.large_run:
                if self.bestqueue not in ("himem.q", "default.q"):
                    ## this queue has no wall clock time limit. make.shared was running out
                    ## of 12 hour limit in medium.q for 3K samples
                    self.bestqueue = "default.q"

### sanity check, this should match the counters                
    def pollrunning(self):
        tmp=defaultdict(int)
        for jid, value in self.registered.values():
            tmp[value]+=1
        self.stats = tmp    
        self.running = self.stats["running"]
    
    def dispatch(self): 
    
        while self.nodesAvailable():
            if not self.scheduled.empty():
                
                tmp = self.scheduled.get()
                self.processing[tmp.getUniqueID()]=tmp
                #print "submitting", tmp.getUniqueID()
                
                
                jid = tmp.submit()
                #print jid
                
                if jid==-1:
                    #whate happens if SGE is temporarily not availabe:
                    #job is marked as running but with
                    #non-existing job ID -1; it is found as "finished" by the
                    #next polling cycle and eventually marked as failed by the
                    #Step thread. However, isJobDone() as it is written will
                    #not mark such -1 job as done as long as SGE is not
                    #available. self.pollqueues() will likely raise in that
                    #case the way it is written now, causing this thread to
                    #terminate.
                    pass
        
                
                self.registered[tmp.getUniqueID()] = [tmp.getGridId(), "running"]   
                self.available-=1
            else:
                break
            
            
        
    def pickQ(self):
        return self.bestqueue
    
    def register(self, task):
    #this is called from Step threads; it could be preemptied
    #by other methods in this thread between the next two lines.
    #The present order should work OK if preemptied by
    #dispatch()
        self.registered[task.getUniqueID()]=[-1, "new"]
        self.scheduled.put(task)
                                
    def shutdown(self):
        self.active=False
        print "Queue status shutting down..."
    
    def nodesAvailable(self):
        return (self.available > 0)
    
    def __str__(self):
        otpt ="Currently running/waiting: %s/%s\n" % (self.running, self.scheduled.qsize())
        otpt ="%savailable/total: %s/%s" % (otpt, self.available, self.maxnodes)
        # for key, tup in self.registered.items():
#            id, status = tup
#            if id != -1:
#                otpt = "%s\n\t%s\t%s\t%s" % (otpt, id, status, key[0:10])
        for key, val in self.stats.items():
            otpt = "%s\n\t%s\t%s" % (otpt, key, val) 

        otpt = "%s\n\nbest queue: %s" % (otpt, self.bestqueue)
        return (otpt)   


    #################################################
    ### a thread that will track of a qsub job
    ### templates adapted to JCVIs grid
    ### 
class GridTask():
    def __init__(self, template="default.q", command = "", name="default", 
            cpu="1", dependson=list(), cwd=".", debug=None, mem_per_cpu=2,
            sleep_start=0, flag_completion=False):

        if debug is None:
            debug = YAPGlobals.debug_grid_tasks

        self.flag_file = None
        
        self.gridjobid=-1
        self.completed=False
        self.queue=template
        self.inputcommand = command
        
        self.cwd=cwd
        self.project = __projectid__
        self.email = __email__
    
        ### remove *e##, *pe## *o## *po##  
        self.retainstreams=" -o /dev/null -e /dev/null "
        
        ### debug flag
        self.debugflag = debug
             
        ncpu = 1
        try:
            ncpu = int(cpu)
        except:
            pass
        
        ### the only queue that has 4 CPUs allowed in pe
        if ncpu>4:
            self.queue = "himem.q"
                
        mem = mem_per_cpu * ncpu

        ## this is the max currently allowed - go figure...
        if mem > 40:
            mem = 40

        if mem > 32:
            self.queue = "himem.q"

        mem_per_cpu = int(mem/ncpu)

        ## -l memory interacts with -pe threaded, resulting in a multiple of them
        ## for the total memory requested
        if mem_per_cpu != 0: 
            mem_spec = "-l memory={}G".format(mem_per_cpu)
        else:
            mem_spec = ""
        
        if len(dependson)>0:
            holdfor = "-hold_jid "
            for k in dependson:
                raise ValueError("Job dependencies do not work yet because jobs are submitted from a separate thread and getGridId() returns -1 at this point")
                holdfor = "%s%s," % (holdfor, k.getGridId())
            holdfor=holdfor.strip(",")      
            sleep_start = max(sleep_start,10) # let files on shared FS to appear
        else:
            holdfor = ""    
        
        ### keep po pe o e streams for debugging purposes
        if self.debugflag:
            self.retainstreams=""
        
            
        ### To source the proper environment, create a script with the command, 
        ### and first source the RC file, then invoke that instead of submitting
        ### the command directly. Note that the -V option to qsub does not
        ### propagate LD_LIBRARY_PATH, which is squashed by the kernel for sudo
        ### programs (which is SGE exec daemon). Note (AT): -V just does not work
        ### for me at all - jobs silently fail (2014-04-10).
        px = "tmp.%s.%s.%s.%s." % (randrange(1,100),randrange(1,100),randrange(1,100),randrange(1,100))
        sx = ".%s.%s.%s.%s.sh" % (randrange(1,100),randrange(1,100),randrange(1,100),randrange(1,100))
        

        ##### to avoid too many opened files OSError
        pool_open_files.acquire()
        ### bounded semaphore should limit throttle the files opening for tasks created around the same time
        try:
            try:
                scriptfile, scriptfilepath = tempfile.mkstemp(suffix=sx, prefix=px, dir=self.cwd, text=True)
            finally:
                os.close(scriptfile)
            self.scriptfilepath = scriptfilepath
            os.chmod(self.scriptfilepath,  0777 )
            if sleep_start > 0:
                import random
                sleep_cmd = "sleep {}".format(sleep_start+random.randint(1,5))
            else:
                sleep_cmd = ""

            input= "#!/bin/bash\n. {}\n{}\n{}\n".format(rcfilepath,sleep_cmd,self.inputcommand)
            
            if flag_completion:
                try:
                    flag_fobj, flag_file = tempfile.mkstemp(suffix="flagfile", prefix=px, dir=self.cwd, text=True)
                    os.write(flag_fobj,"Created\n")
                finally:
                    os.close(flag_fobj)
                self.flag_file = flag_file
                
                flag_file_base = os.path.basename(self.flag_file)
                input += """\nif [[ "$?" != "0" ]]; then (echo {check_flag_str_Fail} > {flag_file_base}); else (echo {check_flag_str_OK} > {flag_file_base}); fi\nsync\n""".\
                        format(flag_file_base=flag_file_base,
                                check_flag_str_OK=check_flag_str_OK,
                                check_flag_str_Fail=check_flag_str_Fail)

            with open(self.scriptfilepath, "w") as scriptfile:
                scriptfile.write(input)
        finally:
            pool_open_files.release()
            ####


        self.templates=dict()
        self.templates["himem.q"]     = 'qsub %s -P %s -N jh.%s -cwd -pe threaded %s %s -l "himem" -M %s -m a  %s "%s" ' % (self.retainstreams, self.project, name, cpu, mem_spec, self.email, holdfor, self.scriptfilepath) 
        self.templates["default.q"]  = 'qsub %s -P %s -N jd.%s -cwd -pe threaded %s %s -M %s -m a %s "%s" ' % (self.retainstreams, self.project, name, cpu, mem_spec, self.email, holdfor,  self.scriptfilepath)
        self.templates["fast.q"]       = 'qsub %s -P %s -N jf.%s -cwd -pe threaded %s %s -l "fast" -M %s -m a  %s "%s" ' % (self.retainstreams, self.project, name,cpu, mem_spec, self.email, holdfor, self.scriptfilepath) 
        self.templates["medium.q"]   = 'qsub %s -P %s -N jm.%s -cwd -pe threaded %s %s -l "medium" -M %s -m a  %s "%s" ' % (self.retainstreams, self.project, name, cpu, mem_spec, self.email, holdfor, self.scriptfilepath)
        self.templates["himemCHEAT"] = 'qsub %s -P %s -N jH.%s -cwd -pe threaded %s %s -l "himem" -M %s -m a  %s "%s" ' % (self.retainstreams, self.project, name, 1, mem_spec, self.email, holdfor, self.scriptfilepath)    
        self.templates["mpi"]        = 'qsub %s -P %s -N jP.%s -cwd -pe orte %s %s -M %s -m a  %s mpirun -np %s "%s" ' % (self.retainstreams, self.project, name, cpu, mem_spec, self.email, holdfor, cpu, self.scriptfilepath )
        self.command = ""
        QS.register(self);
                        
    def submit(self):

        if not self.queue in self.templates.keys():
            self.queue = QS.pickQ()
        self.command = self.templates[self.queue]
        BOH.toPrint("-----","BATCH","command: '{}'".format(self.command))
        
        if YAPGlobals.dummy_grid_tasks:
            check_call(["bash",self.scriptfilepath],cwd=self.cwd, close_fds=True)
            self.gridjobid = 0
        else:
            err = ""
            for i_try in range(3):
                time.sleep(2**(i_try+1)-2)
                p = Popen(shlex.split(self.command), stdout=PIPE, stderr=PIPE, cwd=self.cwd, close_fds=True)
                out, err = p.communicate()
                    
                err = err.strip()
                out = out.strip()
                
                if p.returncode == 0:
                    assert out.endswith("has been submitted"),"Unexpected content in qsub output: {}".format(out)
                    self.gridjobid = out.split(" ")[2]
                    BOH.toPrint("-----","BATCH","qsub output '{}'".format(out))
                    BOH.toPrint("-----","BATCH","getGridId {}".format(self.getGridId()))
                    break
                BOH.toPrint("-----","BATCH","qsub error '{}', trying again...".format(err))
            else:
                BOH.toPrint("-----","BATCH","submit() multiple qsub errors '{}', giving up".format(err))

                    
        return (self.getGridId())

    
    def getGridId(self):
        return self.gridjobid   
    
    def getUniqueID(self):
        return "%s_%s_%s" % (id(self), self.cwd, self.inputcommand)

    def setCompleted(self):
        try:
            if not self.debugflag:
                os.remove(self.scriptfilepath)
        except OSError, error:
            print( "%s already gone" % self.scriptfilepath)
        QS.flagRemoval(self)
        with QS.any_task_completed:
            self.completed = True
            QS.any_task_completed.notify_all()

    def isCompleted(self):
        return self.completed
    
    def wait(self):
        with QS.any_task_completed:
            while not self.isCompleted():
                QS.any_task_completed.wait()
            if self.flag_file:
                return check_flag_file(self.flag_file)
            else:
                return True
 

    #################################################
    ### Iterator over input fasta file.
    ### Only reading when requested
    ### Useful for very large FASTA files
    ### with many sequences
class   FastaParser:
    def __init__ (self, x):
        self.filename = x
        self.fp = open(x, "r")  
        self.currline = "" 
        self.currentFastaName = ""
        self.currentFastaSequence = ""
        self.lastitem=False 
            
    def __iter__(self):
        return(self)    
                
        ##### 
    def next(self):
        for self.currline in self.fp:
            if self.currline.startswith(">"):
                self.currline = self.currline[1:]
                if self.currentFastaName == "":
                    self.currentFastaName = self.currline
                else:
                    otpt = (self.currentFastaName.strip(), self.currentFastaSequence.strip())
                    self.currentFastaName = self.currline
                    self.currentFastaSequence = ""  
                    self.previoustell = self.fp.tell()
                    return (otpt)
                
            else:
                self.addSequence(self.currline) 
        
        if not self.lastitem:
            self.lastitem=True          
            return (self.currentFastaName.strip(), self.currentFastaSequence.strip())
        else:
            raise StopIteration 
                                
    def addSequence(self, x):
            self.currentFastaSequence = "%s%s" % (self.currentFastaSequence, x.strip())         
                        
    def __str__():
        return ("reading file: %s" %self.filename)  

class GeneralPurposeParser:
    def __init__(self, file, skip=0, sep="\t", skip_empty=True):
        self.skip_empty = skip_empty
        self.filename = file
        self.fp = open(self.filename, "rU") 
        self.sep = sep
        self.skip = skip
        self.linecounter = 0
        self.currline=""
        while self.skip>0:
            self.next()
            self.skip-=1
        
    def __iter__(self):
        return (self)
    
    def next(self):
        for currline in self.fp:
            currline = currline.strip()
            self.linecounter = self.linecounter + 1
            if not self.skip_empty or currline:
                currline = currline.split(self.sep)
                self.currline = currline
                return(currline)            
        raise StopIteration
                    
    def __str__(self):
        return "%s [%s]\n\t%s" % (self.filename, self.linecounter, self.currline)

    #################################################
    ### The mother of all Steps:
    ###

if YAPGlobals.step_dummy_thread:
    DefaultStepBase = ReportingDummyThread
else:
    DefaultStepBase = ReportingThread

class   DefaultStep(DefaultStepBase):
    #This limits the number of concurrent ("submitted") step instances,
    #in other words, the number of threads on which start() has been 
    #called. The code creating instances of derived classes will block
    #until the semaphore is acquired.
    semaphore = BoundedSemaphore(YAPGlobals.step_threads_max)
    def __init__(self):
        #### thread init
        DefaultStepBase.__init__(self)
        self.random = uniform(0, 10000)
        self.name = ("%s[%s]" % (self.name, self.random))

        #### hash of the current step-path (hash digest of previous steps + current inputs + arguments + control)
        self.workpathid = None
        #### path where the step stores its files
        self.stepdir = ""
                    
        #### what needs to be completed for this step to proceed
        #### a list of steps
        self.previous = list()
        
        #### mapping type - path for files
        self.inputs = defaultdict(set)
        
        #### mapping type - name for files
        self.outputs = defaultdict(set)
        
        #### mapping arg  val for program's arguments
        self.arguments= dict()
        
        #### mapping arg  val for step's control parameters
        self.control= dict()
        
        #### ID of the step...
        self.stepname = ""  
        
        #### flag for completion
        self.completed = False
        self.completedpreviously=False
        self.failed = False
        
        #### keep track of time elapsed
        self.starttime = 0
        self.endtime = 0
        
        #### special flag, some steps might not want to delete the inputs (argcheck)
        self.removeinputs = True
        
        ####
        
    def setInputs(self, x):
        for k,v in x.items():
            for elem in v:
                self.inputs[k].add(elem)
    
    def setArguments(self, x):
        for k,v in x.items():
            if v=="":
                v=" "
            self.arguments[k] = v
    
    def setControl(self, x):
        for k,v in x.items():
            if v=="":
                v=" "
            self.control[k] = v

    def setPrevious(self, x):
        if isinstance(x,DefaultStep):
            self.previous.append(x)
        else:
            for elem in x:
                
                if not isinstance(elem,DefaultStep):
                    msg = "List of PREVIOUS steps in step {} must contain Steps objects, instead found this element: {}".\
                            format(self.getName(),elem)
                    self.message(msg)
                    raise ValueError(msg)

                self.previous.append(elem)
    
    def setName(self, x):
        self.stepname=x

    def getName(self):
        return self.stepname
    
    #override the start method to make sure that we do not
    #generate any exceptions in our code between acquiring semaphore
    #and calling Thread.start()
    def start(self):
        DefaultStep.semaphore.acquire()
        try:
            DefaultStepBase.start(self)
        except:
            DefaultStep.semaphore.release()
            raise
    
    def do_run(self):
        try:
            try:
                self.init()
            except:
                traceback.print_exc()
                self.message(traceback.format_exc())
                self.failed = True

            if self.failed:
                ##deregister never throws
                BOH.deregister()
                self.completed=True
            elif not self.isDone():
                try:
                    self.performStep()
                    self.finalize()
                except:
                    traceback.print_exc()
                    self.message(traceback.format_exc())
                    ##deregister never throws
                    BOH.deregister()
                    self.completed=True
                    self.failed=True
            else:
                self.message("Completed (previously).") 
                BOH.deregister()
                self.completed=True
                self.completedpreviously=True
        finally:
            DefaultStep.semaphore.release()
        
    def performStep(self):
        self.message("in a step...")
                
    def init(self):
    
        redo=False
        ### wait for previous steps to finish
        for k in self.previous:
            k.join()
            #while not k.isDone():
            #    time.sleep(0.01)
            if k.hasFailed():
                self.failed=True
            redo=redo or (not k.isDonePreviously()) 
        
        #self.message("needs a redo %s" % (redo))
        if not self.failed:
            ### time stamp
            self.starttime = time.time()
        
            #### hash of the current step-path (hash digest of previous steps + current inputs + arguments?)
            self.workpathid = self.makeWorkPathId()
            ####
            
            ### output handler
            BOH.register(self.workpathid)
            ###
            
            #self.message("Initializing %s %s" % (self.workpathid, self.name))
            
            #### create directories if necessary
            self.stepdir =""
            self.prepareDir(redo=redo)  
            
    def makeWorkPathId(self):
        tmp = list()
        tmp.append(self.stepname)
        if self.previous!=None:
            for k in self.previous:
                while not k.getWorkPathId():
                    time.wait(0.1)
                tmp.extend([k.getWorkPathId()])
            
        for k,v in self.inputs.items():
            tmp.extend(["%s=%s" % (k, ",".join(v) ) ] )
        
        for k,v in self.arguments.items():
            tmp.extend(["%s=%s" % (k, v) ] )
                
        for k,v in self.control.items():
            tmp.extend(["%s=%s" % (k, v) ] )
                
        tmp.sort()  
        tmp = "\n".join(tmp)
        
        workpathid = hashlib.md5(tmp).hexdigest()
        return (workpathid)
                
    def getWorkPathId(self):    
        return (self.workpathid)
                
    def prepareDir(self, redo=False):
        assert self.workpathid
        ### make step's directory
        self.stepdir = "Step_%s_%s" % (self.stepname, self.workpathid)
        
        
        flush_old = False
        try:    
            os.mkdir(self.stepdir)              
        except OSError, error:
            self.message( "Step directory already exists...")
            flush_old=True
        
        
        if redo:
            if flush_old:
                self.message("Updating...")
                k = "rm -r *"
                task = GridTask(template="pick", name="redo_clean", command=k, cpu=1,  cwd = self.stepdir)
                task.wait()
            else:
                ###supposedly no old step-data to flush
                pass    
            
        else:   
            ### has analysis been done already?
            try:
                if self.parseManifest():
                    self.completed=True
                    self.completedpreviously=True
                    self.message("Using data generated previously...")
            except: 
                traceback.print_exc()
                self.message("****ERROR***")
                self.message(traceback.format_exc())
                self.message("************")                    
                                                                    
    def finalize(self):
        
        if not self.failed:         
            self.categorizeAndTagOutputs()
            self.makeManifest()
            
            self.endtime = time.time()
            self.message( "+%s\t[Done]" % (str(datetime.timedelta(seconds=round(self.endtime-self.starttime,0))).rjust(17)) ) 
        else:
            self.endtime = time.time()
            self.message( "+%s\t[Fail]" % (str(datetime.timedelta(seconds=round(self.endtime-self.starttime,0))).rjust(17)) ) 
            
        self.completed=True 
        BOH.deregister()

    def waitOnGridTasks(self,tasks,fail_policy="ExitOnFirst"):
        fail_policy_values = ('ExitOnFirst','ExitOnLast','ExitNever','Ignore')
        assert fail_policy in fail_policy_values, "Allowed values: {}".format(fail_policy_values)
        failed_any = False
        failed_descr = ""
        for task in tasks:
            res = task.wait()
            failed_any |= (not res)
            if fail_policy != "Ignore" and not res:
                self.fail = True
                failed_descr = "Grid ID: {}. Descr: {}".format(task.getGridId(),getUniqueID())
                if fail_policy == "ExitOnFirst":
                    raise ValueError("Step {} had grid task failing, aborting step. {}".format(self.stepname,failed_descr))
        
        if fail_policy == "ExitOnLast" and failed_any:
            raise ValueError("Step {} had at least one grid task failing, aborting step. {}".format(self.stepname,failed_descr))
        return not failed_any

    def makeManifest(self):
        assert self.workpathid
        pool_open_files.acquire()
        try:
            with open("%s/%s.manifest" % \
                  (self.stepdir, self.workpathid), "w") as m:
                for type, files in self.inputs.items():
                    if len(files)>0:
                        m.write("input\t%s\t%s\n" % (type, ",".join(files)) )
                for arg, val in self.arguments.items():
                    m.write("argument\t%s\t%s\n" % (arg, val ) )
                for arg, val in self.control.items():
                    m.write("control\t%s\t%s\n" % (arg, val ) )
                for type, files in self.outputs.items():
                    if len(files)>0:
                        m.write("output\t%s\t%s\n" % (type, ",".join(files)) )
        finally:
            pool_open_files.release()
    
    def determineType(self, filename):
        filename = filename.strip().split(".")
        extension = filename[-1]
        preextension = filename[-2]
        
        if "scrap" in filename and extension == "fasta": #scrap.contigs.fasta comes from make.contigs
            return "scrapfasta"
        
        elif "scrap" in filename: #scrap.contigs.fasta comes from make.contigs
            return "scrap"
        
        elif extension == "refalign":
            return "refalign"
            
        elif preextension == "align" and extension == "report":
            return "alignreport"
    
        elif extension == "dist" and preextension == "phylip":
            return "phylip"         
        elif extension == "dist":
            return "column"
                
        elif preextension == "tax" and extension =="summary":
            return "taxsummary"
        
        elif preextension == "cdhit" and extension =="clstr":
            return "cdhitclstr"
        elif preextension == "bak" and extension =="clstr":
            return "cdhitbak"
        elif extension == "cdhit":
            return "fasta"
        
        elif extension in ["align", "fna", "fa", "seq", "aln"]:
            return "fasta"
            
        elif extension == "qual":
            return "qfile"
        
        elif extension == "tax":
            return "taxonomy"
        
        elif extension in ("count_table","count"):
            return "count"
        elif extension == "names":
            return "name"
        elif extension == "groups":
            return "group"
        elif extension == "files":
            return "file"
        
        elif extension in ["tre", "tree", "dnd"]:
            return "tre"
        
        ### sge job files
        elif re.search("po\d{3}", extension) != None:
            return "po"
        elif re.search("pe\d{3}", extension) != None:
            return "pe"     
        elif re.search("o\d{3}", extension) != None:
            return "o"
        elif re.search("e\d{3}", extension) != None:
            return "e"  
        else:
            return extension
        
    def categorizeAndTagOutputs(self):
        assert self.workpathid
        inputs = [x.split("/")[-1] for x in unlist( self.inputs.values()) ]
        for file in glob.glob("%s/*" % self.stepdir):
            file = file.split("/")[-1]
            if file in inputs:
                if self.removeinputs:
                    command = "unlink %s" % (file)
                    p = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, cwd=self.stepdir, close_fds=True) 
                    out,err = p.communicate()
                    

                else:
                    self.message("kept %s" % file)
                    #pass
            elif not file.endswith("manifest"):
                #self.message( "output: %s" % (file))
                
                ### make sure that every output file except for the manifest starts with the workpathID
                file = file.split(".")
                if len(file[0]) == len(self.workpathid):
                    newfilename = "%s.%s"  % (self.workpathid, ".".join(file[1:]))
                else:
                    newfilename = "%s.%s"  % (self.workpathid, ".".join(file[0:])) 
                
                
                if ".".join(file) != newfilename:
                    k="mv %s %s" % (".".join(file), newfilename)
                    p = Popen(shlex.split(k), stdout=PIPE, stderr=PIPE, cwd=self.stepdir, close_fds=True)   
                    out,err = p.communicate()
                    

                self.outputs[self.determineType(newfilename)].add(newfilename)
                        
    def find(self, arg, ln=True, original=False, require=None):
        if require is None:
            require = YAPGlobals.find_require
        files=list()   
        if not original:        
            if len(self.inputs[arg])==0:
                tmp = {arg: self.getOutputs(arg)}
                self.setInputs(tmp)     
        else:
            tmp = {arg: self.getOriginal(arg)}      
            self.setInputs(tmp) 

        files = self.inputs[arg]
        if require and not len(files):
            self.message("No files found that match {} and 'require' is set. Raising an exception.".format(arg))
            raise ValueError(arg)
        
        toreturn=list() 
        
        for file in files:
            if self.isVar(file):
                toreturn.append(file[5:])
            else:   
                targ = os.path.basename(file)
                targ_path = os.path.join(self.stepdir,targ)
                rmrf(targ_path)
                if (ln):
                    command = "cp -s %s %s" % (file, targ )
                else:
                    command = "cp %s %s" % (file, targ )
                check_call(shlex.split(command), cwd=self.stepdir, close_fds=True )    
                toreturn.append(targ)        
        #unique 
        toreturn = set(toreturn)
        return list(toreturn)
                
    def isVar(self,x):
        return x.startswith("[var]")
    

    def setOutputMaskInclude(self,patt):
        """Set output inclusion pattern for lookups by upstream steps.
        patt must be Python regex. type1|type2|type3 will match either one of three types.
        """
        self.control["output_mask_include"] = patt
    
    def setOutputMaskExclude(self,patt):
        """Set output exclusion pattern for lookups by upstream steps.
        patt must be Python regex. type1|type2|type3 will match either one of three types.
        """
        self.control["output_mask_exclude"] = patt

    def isMasked(self,key):
        """Return True if control parameters mask this data type from lookups by downstream steps.
        include pattern has preference over exclude pattern"""
        incl = self.control.get("output_mask_include",None)
        if incl is not None:
            if re.match("^"+incl+"$",key):
                return False
        excl = self.control.get("output_mask_exclude",None)
        if excl is not None:
            if re.match("^"+excl+"$",key):
                return True
        return False

    
    def getOutputs(self, arg):
        
        otpt = list()
        
        is_masked = self.isMasked(arg)

        if not (self.isDone() and is_masked):
            ## when step itself is executing, it is allowed
            ## to ask upstream for types that it masks
        
            if self.outputs.has_key(arg):
                for x in unlist(self.outputs[arg]):
                    if self.isVar(x):
                        otpt.append(x)
                    else:
                        otpt.append("../%s/%s" % (self.stepdir, x))
                
            elif self.previous!=None:
                    for k in self.previous:
                        otpt.extend(k.getOutputs(arg))
        
        return otpt
        
    def getOriginal(self, arg):
        if self.previous == None:
            return self.getOutputs(arg)
        else:   
            current = self.getOutputs(arg)
            otpt = list()           
            for k in self.previous:
                otpt.extend(k.getOriginal(arg))
            if len(otpt)>0:
                return otpt
            else:
                return current      
            
    def parseManifest(self):
        assert self.workpathid
        manifest_file = "%s/%s.manifest" % (self.stepdir, self.workpathid)
        if not os.path.exists(manifest_file):
            return False
        pool_open_files.acquire()
        try:
            with open(manifest_file,"r") as fp:
                lines=fp.readlines()
        finally:
            pool_open_files.release()
        for line in lines:
            line = line.strip("\n").split("\t")
            if line[0] == "output":
                type = line[1]
                files = line[2].split(",")
                for file in files:
                    self.outputs[type].add(file)
            elif line[0] == "input":
                type = line[1] 
                files = line[2].split(",")
                for file in files:
                    self.inputs[type].add(file) 
            elif line[0] == "argument":
                if len(line)==2:
                    self.arguments[line[1]] = " "
                else:
                    self.arguments[line[1]]=line[2]
            elif line[0] == "control":
                if len(line)==2:
                    self.control[line[1]] = " "
                else:
                    self.control[line[1]]=line[2]
        return True
                
    def message(self, text):
        date = time.localtime()
        if type(text) == list:
            for line in text: 
                BOH.toPrint(self.workpathid, self.stepname, line, date)
        else:   
            BOH.toPrint(self.workpathid, self.stepname, text, date)       
                
    def isDone(self):
        return self.completed 
    
    def isDonePreviously(self):
        return self.completedpreviously 

    def hasFailed(self):
        return self.failed
        
    def getInputValue(self, arg):
        if self.arguments.has_key(arg):
            return self.arguments[arg]
        else:
            return None
            
    def setOutputValue(self, arg, val):
        self.outputs[arg] = ["[var]%s" % (val)] 
                
                
    def __str__(self):
        otpt = "%s\t%s" % (self.stepname, self.name)
            
        for val in self.previous:
            otpt += "%s\n%s" % (otpt, val.__str__())    
            
        #otpt = "\n".join(set(otpt.strip().split("\n")))        
        return otpt 
    
class   FileImport(DefaultStep):    
    def __init__(self, INS):
        DefaultStep.__init__(self)
        self.setInputs(INS)
        #self.setArguments(ARGS)
        #self.setPrevious(PREV)
        self.setName("FILE_input")
        self.start()
            
    def performStep(self):
        for type in self.inputs.keys():
            files = self.inputs[type]
            for file in files:
                pool_open_files.acquire()
                try:
                    file = file.split("~")
                    if len(file)>1:
                        file, newname = file
                        tmp = file.strip().split("/")[-1]
                        k = "cp %s %s.%s" % (file, newname, type)
                    else:
                        file = file[0]
                        k ="cp %s imported.%s.%s"  % (file, os.path.basename(file), type)
                    self.message(k)
                    subprocess.call(shlex.split(k), cwd=self.stepdir, close_fds=True)
                finally:
                    pool_open_files.release()
                
class   ArgumentCheck(DefaultStep):
    def __init__(self, SHOW, PREV):
        ARGS = {"show":SHOW}
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("ArgCheck")
        #self.nodeCPUs=nodeCPUs
        self.removeinputs=False
        self.start()
        
    def performStep(self):
        x = self.getInputValue("show")
        if x!=None:
            for type in x.split(","):
                for file in self.find(type):
                    self.message("%s: %s" % (type,file))    

class   OutputStep(DefaultStep):
    def __init__(self, NAME, SHOW, PREV):
        ARGS = {"show":SHOW}
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("OUTPUT_%s" % (NAME))
        #self.nodeCPUs=nodeCPUs
        self.removeinputs=False
        self.start()
        
    def performStep(self):
        x = self.getInputValue("show")
        if x!=None:
            for type in x.split(","):
                for file in self.find(type.strip(), ln = False):
                    self.message("%s: %s" % (type,file))
        
class   SFFInfoStep(DefaultStep):
    def __init__(self, INS, ARGS, PREV):
        DefaultStep.__init__(self)
        self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("sffinfo")
        self.start()
        
    def performStep(self):
        steps = list()      
        for sff in self.find("sff"):
            
            k = "/usr/local/bin/sffinfo -s %s > %s.fasta" % (sff, sff)
            self.message(k)
            task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir)
            steps.append(task)
            
            k = "/usr/local/bin/sffinfo -q %s > %s.qual" % (sff, sff)
            self.message(k)
            task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir)
            steps.append(task)
            
            k = "/usr/local/bin/sffinfo -f %s > %s.flow" % (sff, sff)
            self.message(k)
            task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir)
            steps.append(task)
            
        for s in steps:
            s.wait()

class   WriteOligosStep(DefaultStep):
    
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("WriteOligos")
        self.start()
        
    def performStep(self):
        out_name = os.path.join(self.stepdir,"seq.oligos")
        with open(out_name,"w") as out:
            out.write("primer {} {}\n".format(
                self.arguments["forward_primer"],
                self.arguments["reverse_primer"])
                )
    
class   WriteFilesStep(DefaultStep):
    """Dump values of arguments dict with file names provided by the keys.
    The use case is reasonably small values, such as writing accnos file
    for a fixed set of sequence names"""
    
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("WriteFiles")
        self.start()
        
    def performStep(self):
        for file_name, file_content in self.arguments.items():
            out_name = os.path.join(self.stepdir,file_name)
            with open(out_name,"w") as out:
                out.write(file_content)
    
class   AlignmentReportTemplateRangeStep(DefaultStep):  
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("AlignmentReportTemplateRange")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
                
    def performStep(self):
        import csv
        template_name = self.arguments.get("template_name",None)
        f = self.find("alignreport",require=True)[0]
        start = 10**10
        end = -1
        template_length = 0
        template_found = False
        with open(f,"r") as inp_file:
            inp = csv.DictReader(inp_file,delimiter="\t")
            for rec in inp:
                if template_name is None:
                    template_name = rec["TemplateName"]
                if template_name == rec["TemplateName"]:
                    start = min(start,int(rec["TemplateStart"]))
                    end = max(end,int(rec["TemplateEnd"]))
                    template_length = rec["TemplateLength"]
                    template_found = True

        assert template_found
        pad = self.arguments.get("template_range_pad",0)
        self.setOutputValue("template_range_start",max(1,start-pad))
        self.setOutputValue("template_range_end",min(template_length,end+pad))

class   PcrReferenceAlignmentStep(DefaultStep):
    def __init__(self, ARGS, PREV):      
        DefaultStep.__init__(self)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("PcrReferenceAlignment")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        call_args = dict(globals())
        call_args.update(self.arguments)
        call_args["reference_ali"] = self.find(self.arguments["align_type"],require=True)[0]
        call_args["pcr_target_fasta"] = self.find("fasta",require=True)[0]
        call_args["mothurpath"] = mothurpath

        k = """{binpath}python {scriptspath}/scripts.py pcr-reference-alignment \
                --pcr-target-padding {pcr_target_padding} \
                {mothurpath}mothur \
                {pcr_target_fasta} \
                {pcr_target_idseq} \
                {primer_forward} \
                {primer_reverse} \
                {reference_ali}""".format(**call_args)
                        
        self.message(k) 
        task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
        task.wait() 


class   MothurStep(DefaultStep):
    def __init__(self, NM, nodeCPUs, INS, ARGS, PREV):
        DefaultStep.__init__(self)
        self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName(NM)
        self.nodeCPUs=nodeCPUs
        self.start()
    
    def makeCall(self):
        skip_step = False
        FORCE = self.getInputValue("force")
        x = MOTHUR.getCommandInfo(self.stepname)
        #self.message(self.inputs)
        if FORCE != None :
            TYPES = [ t.strip() for t in FORCE.strip().split(",") if t.strip() != "" ]
        else:
            TYPES = x.getInputs()
        
        FORCE_EXCL = self.getInputValue("force_exclude")
        if FORCE_EXCL != None:
            TYPES_EXCL = FORCE_EXCL.strip().split(",")
            TYPES = [ t for t in TYPES if t not in TYPES_EXCL ]
            
        mothurargs = list()
        
        for arg, val in self.arguments.items():

            if arg not in ("force","force_exclude"):
        
                if x.isAnArgument(arg):
                    mothurargs.append("%s=%s" % (arg, val))
                elif arg in ("find","find_req","find_or_skip_step"):
                    for a in val.strip().split(","):
                        self.message(a)
                        kv = a.strip().split("=",1)
                        if len(kv) > 1:
                            k,v = kv
                        else:
                            k = kv[0]
                            v = k
                        v = v.strip()
                        k = k.strip()
                        valstoinsert = self.find(v)
                        self.message(valstoinsert)
                        if valstoinsert:
                            TYPES = [ t for t in TYPES if t != k ]
                        valstoinsert = [ val_el.strip() for val_el in valstoinsert if val_el.strip() ]
                        if valstoinsert:
                            mothurargs.append("%s=%s" % (k, "-".join(valstoinsert)))
                        else:
                            if arg == "find_req":
                                self.message("Required argument '{}' not found!".format(v))
                                raise ValueError(v)
                            elif arg == "find_or_skip_step":
                                self.message("Argument '{}' was marked for skipping Mothur step if not found; skipping step.".format(v))
                                skip_step = True
                            else:
                                self.message("Skipping Mothur argument insertion '{}' - not found".format(v))
                else:
                    self.message("skipping '%s', as it is not an argument for %s" % (arg, self.stepname))
        
        for TYPE in TYPES:
            #### on occasion, mothur overwrites the original file - namely names file
            #### FALSE creates a copy
            #### TRUE creates a link
            if TYPE=="name":
                tmp = self.find(TYPE, False)
            else:
                tmp = self.find(TYPE, True)
                
            if len(tmp)>0:
                mothurargs.append ("%s=%s" % (TYPE, "-".join(tmp)))
            else:
                if x.isRequired(TYPE):
                    self.message("Required argument '%s' not found!" % (TYPE))  
                    raise ValueError(TYPE)
                else:
                    self.message("Optional argument '%s' not found, skipping" % (TYPE)) 
        
    
        ### method is parallelizable,   
        if x.isAnArgument("processors") and "processors" not in self.arguments.keys():
            mothurargs.append("%s=%s" % ("processors", self.nodeCPUs ))
            self.message("Will run on %s processors" % (self.nodeCPUs))
        
        himemflag=False
        ### steps requiring lots of memory
        if self.stepname in ("clearcut", "align.seq"):
            himemflag=True
            self.message("Needs lots of memory")

        if skip_step:
            command = None
        else:
            command = "%s(%s)" % (self.stepname, ", ".join(mothurargs))
        
        return (command, x.isAnArgument("processors"), himemflag)
    
    def performStep(self):
        
        call, parallel, himem = self.makeCall() 

        if call:

            k = "%smothur \"#%s\"" % (mothurpath, call)
            
            self.message(k)
            if (parallel and self.nodeCPUs>1):
                #task = GridTask(template=defaulttemplate, name=self.stepname, command=k, cpu=self.nodeCPUs, dependson=list(), cwd = self.stepdir)
                task = GridTask(template="pick", name=self.stepname, command=k, cpu=self.nodeCPUs, dependson=list(), cwd = self.stepdir)
            
            #elif (himem):
            #   task = GridTask(template="himem.q", name=self.stepname, command=k, cpu=self.nodeCPUs, dependson=list(), cwd = self.stepdir)
            else:
                task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir)
            
            task.wait()
            self.parseLogfile()
        
        else:

            self.message("Step was marked to be skipped, skipping")
        
    def parseLogfile(self):
        for f in glob.glob("%s/*.logfile" % (self.stepdir)):
            line = ""
            for line in loadLines(f):
                ### UCHIME throws an error when it does not find chimeras, even though it completes.
                if line.find ("ERROR")>-1 and line.find("uchime")==-1:
                    self.failed=True
            
            ### last line
            if line.find("quit()")==-1:
                self.failed=True

class   MothurSHHH(DefaultStep):
    def __init__(self,  PREV, nodeCPUs):        
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        #self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("MPyro")
        ####
        ####self.nodeCPUs=nodeCPUs
        self.nodeCPUs=4
        self.start()
        
    def performStep(self):
        tasks = list()
        
        TOC = self.find("file")
        flows = self.find("flow")
        
        TOC = loadLines("%s/%s" % (self.stepdir, TOC[0]))
        TOC = [ ".".join(x.strip().split(".")[1:]) for x in TOC]
            
#       for f in flows:
#           tmp = ".".join(f.split(".")[1:])
#           
#           if tmp in TOC:
#               
#               ### split tmp into 10,000 lines chunks
#               k = "split -l 7000 -a 3 %s %s.split." % (f, f)
#               task = GridTask(template="pick", name="MPyroSplit", command=k, cpu=1,  cwd = self.stepdir, debug=False)
#               tasks.append(task)              
#           else:
#               self.message("skipping %s" % (f))  
#       
#       self.message("splitting %s file(s)" % len(tasks))
#       
#       for task in tasks:
#           task.wait() 
        
        ################################################
        tasks = list()
        
        #for chunk in glob.glob("%s/*.split.*" % (self.stepdir)):
        #   chunk = chunk.split("/")[-1]
            #self.message(chunk)    
        #   call = "shhh.flows(flow=%s, processors=%s, maxiter=100, large=10000)" % (chunk, self.nodeCPUs)
        for f in flows:
            tmp = ".".join(f.split(".")[1:])
            if tmp in TOC:
                call = "shhh.flows(flow=%s, processors=%s, maxiter=100, large=10000)" % (f, self.nodeCPUs)  
                k = "%smothur \"#%s\"" % (mothurpath, call)
                self.message(k)
                task = GridTask(template="pick", name="Mpyro", command=k, cpu=self.nodeCPUs,  cwd = self.stepdir, debug=True)
                tasks.append(task)
        if len(tasks)==0:
            self.failed=True        
        self.message("processing %s file(s)" % len(tasks))
        
        for task in tasks:
            task.wait()
                            
class   LUCYcheck(DefaultStep):
    def __init__(self, nodeCPUs, PREV):
        DefaultStep.__init__(self)
        self.nodeCPUs=nodeCPUs
        #self.setInputs(INS)
        #self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("LUCY_check")
        self.nodeCPUs=nodeCPUs
        if self.nodeCPUs>32:
            self.nodeCPUs=30
        self.start()
                
    def performStep(self):
        f = self.find("fasta")[0]
        q = self.find("qfile")[0]
        
        statinfo = os.stat("%s/%s" % (self.stepdir, f))
        #self.message(statinfo.st_size)
        
        if statinfo.st_size==0:
            self.message("%s is empty." % f)   
            self.failed=True
        else:
            
            k ="%s/lucy -error 0.002 0.002 -bracket 20 0.002 -debug -xtra %s -output %s.fastalucy %s.qfilelucy %s %s" % (binpath, self.nodeCPUs, f,q, f,q)
            self.message(k)
            if self.nodeCPUs>2:
                task = GridTask(template=defaulttemplate, name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir)
            else:
                task = GridTask(template="pick", name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir)
            task.wait() 
        
class   LUCYtrim(DefaultStep):  
    def __init__(self, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        #self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("LUCY_trim")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
                
    def performStep(self):
        f = self.find("fastalucy")[0]
        q = self.find("qfilelucy")[0]
        
        k = "%spython %s/fastAtrimmer.py -l %s %s %s " % (binpath, scriptspath, f.split(".")[0], f, q)
        self.message(k)
        task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
        task.wait()

class   MatchGroupsToFasta(DefaultStep):
    def __init__(self, INS, PREV):      
        DefaultStep.__init__(self)
        self.setInputs(INS)
        #self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("MatchGroups")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        tasks = list()
        f = self.find("fasta")
        f = f[0]
        g = self.find("group")
        g = g[0]
        
        n = self.find("name")
        if len(n)>0:
            n = "-n %s" % (n[0])
        else:
            n = ""
        
        k = "%spython %s/MatchGroupsToFasta.py %s -f %s -g %s -o %s.matched.group" % (binpath, scriptspath, n, f, g, ".".join(g.split(".")[:-1]))
        self.message(k) 
        task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
        task.wait() 

class   MatchGroupsToList(DefaultStep):
    def __init__(self, INS, PREV):      
        DefaultStep.__init__(self)
        self.setInputs(INS)
        #self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("MatchGroups")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        tasks = list()
        f = self.find("list")
        f = f[0]
        g = self.find("group")
        g = g[0]

        
        k = "%spython %s/MatchGroupsToFasta.py -l %s -g %s -o %s.matched.group" % (binpath, scriptspath, f, g, ".".join(g.split(".")[:-1]))
        self.message(k) 
        task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
        task.wait()         

def sort_strings_by_regex_list(sl,rl):
    """Example:
    sl = ["aaa.0.1.shared","aaa.0.01.shared","aaa.0.03.shared"]
    rl = [r"0.01",r"0.03",r"0.1"]
    print sort_strings_by_regex_list(sl,[re.escape(r) for r in rl])
    >>> ["aaa.0.01.shared","aaa.0.03.shared","aaa.0.1.shared"]
    """
    import re
    
    def rx_key(rl,s):
        for (r_ind,r) in enumerate(rl):
            if re.search(r,s):
                return r_ind
        return len(rl)
    
    return sorted(sl,key = lambda s: rx_key(rl,s))

class   FileMerger(DefaultStep):
    def __init__(self, TYPES, PREV, prefix="files", cut_header_lines_others=0, order=""):
        """@param order is a list of regex patterns; currently ignored when > 25 files w/o header cutting"""
        ARGS =  dict(types=TYPES,
                prefix=prefix, 
                cut_header_lines_others=cut_header_lines_others, 
                order=order) 
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("FILE_cat")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        
        prefix=self.getInputValue("prefix") 
        cut_header_lines_others=self.getInputValue("cut_header_lines_others") 
        order=self.getInputValue("order")
        if order == "":
            order = None

        tasks = list()
        for t in self.getInputValue("types").strip().split(","):
            files = self.find(t)
            #self.message("FileMerger: files found {}".format(files))
            if order is not None:
                files = sort_strings_by_regex_list(files,order)
            #self.message("FileMerger: files after ordering {}".format(files))

            if cut_header_lines_others == 0:
                if len(files)>0 and len(files)<25:
                    k = "cat %s > %s.x%s.merged.%s" % (" ".join(files), prefix, len(files), t)
                    self.message(k) 
                    task = GridTask(template="pick", name="cat", command=k, cpu=1,  cwd = self.stepdir)
                    tasks.append(task)
                elif len(files)>=25:
                    k = "cat *.%s* > %s.x%s.merged.%s" % (t, prefix, len(files), t)
                    self.message(k) 
                    task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
                    tasks.append(task)
            else:
                if len(files) > 0:
                    fn_out = "%s.x%s.merged.%s" % (prefix, len(files), t)
                    k = "cat %s > %s\n" % (files[0], fn_out)
                    ## we want to maintain order of concatenation regarless of the number of files,
                    ## hence the loop instead of the glob
                    for fn_inp in files[1:]:
                        k += "tail -q -n +{} {} >> {}\n".format(cut_header_lines_others+1,fn_inp,fn_out)
                    self.message(k) 
                    task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
                    tasks.append(task)
            #else:
            #   self.failed=True
            
        for task in tasks:
            task.wait()
            
class   FileSort(DefaultStep):
    def __init__(self, TYPES, PREV):
        ARGS =  {"types": TYPES}        
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("FILE_sort")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        tasks = list()
        for t in self.getInputValue("types").strip().split(","):
            files = self.find(t)
            if len(files)>0:
                k = "sort -n %s > files_x%s.sorted.%s" % (" ".join(files), len(files), t)
                self.message(k) 
                task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
                tasks.append(task)
        for task in tasks:
            task.wait()
                        
class   FileType(DefaultStep):
    def __init__(self, ARGS, PREV):     
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("FILE_type")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        tasks = list()
        for input, output in self.arguments.items():
            files = self.find(input)
            for file in files:
                outname = "%s.%s" % (file, output)
                k = "cp %s %s" % (file, outname)
                self.message(k) 
                task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
                tasks.append(task)
        for task in tasks:
            task.wait()

class FileTypeMaskInput(FileType):
    """Type renamer class that will also mask its input types from dependant tasks.
    This is an equivalent of 'move' command in this immutable workflow pattern."""
    ##TODO: currently imposible to call the __init__ of parent Step classes because they
    ##call threading.start() at the end of __init__(). Factor out the actual init into a do_init()
    ##method that should be overaloaded and chain-called in children classes
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("FILE_type_mask_input")
        excl = "|".join( (input for (input,output) in ARGS.items()) )
        self.setOutputMaskExclude(excl)
        self.start()


class   MaskType(DefaultStep):
    """This masks a list of file types so that they are never found by downstream steps"""

    def __init__(self, TYPES, PREV):     
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        excl = "|".join( TYPES.strip().split(",") )
        self.setOutputMaskExclude(excl)
        self.setPrevious(PREV)
        self.setName("MASK_type")
        #self.nodeCPUs=nodeCPUs
        self.start()  
            
class   MaskExceptType(DefaultStep):
    """This masks all file types but the given ones so that they are never found by downstream steps"""

    def __init__(self, TYPES, PREV):     
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        incl = "|".join( TYPES.strip().split(",") )
        self.setOutputMaskInclude(incl)
        self.setOutputMaskExclude(".*")
        self.setPrevious(PREV)
        self.setName("MASK_except_type")
        #self.nodeCPUs=nodeCPUs
        self.start()  
            
class   CleanFasta(DefaultStep):
    def __init__(self, INS, PREV):      
        DefaultStep.__init__(self)
        self.setInputs(INS)
        #self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("CleanFasta")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        tasks = list()
        f = self.find("fasta")
        f = f[0]
        k = "%spython %s/CleanFasta.py -i %s -o %s.dash_stripped.fasta" % (binpath, scriptspath,f, ".".join(f.split(".")[:-1]))
        self.message(k) 
        task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
        task.wait()         
        
class   MakeNamesFile(DefaultStep):
    def __init__(self, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        #self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("FILE_names")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        files = self.find("fasta")
        for f in files:
            self.message("Creating 'names' file for sequences in {0}".format( f))
    
            newname = f.strip().split(".")[:-1]
            newname = "%s.names" % (".".join(newname)) 
            otpt = open("%s/%s" % (self.stepdir,newname ), 'w')
            for head, seq in FastaParser("%s/%s" % (self.stepdir, f)):
                head = head.strip().split()[0]
                otpt.write("%s\t%s\n" % (head, head))
            otpt.close()    
            
        if len(files)==0:
            self.message("No files to generate NAMES...")       
            
class   MakeGroupsFile(DefaultStep):
    def __init__(self, PREV, id):
        ARGS =  {"groupid": id}
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("FILE_groups")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        files = self.find("fasta")
        for f in files:
            id = self.getInputValue("groupid")
            self.message("Creating 'groups' file; '{0}' for sequences in {1}".format(id, f))
            newname = f.strip().split(".")[:-1]
            newname = "%s.groups" % (".".join(newname)) 
            otpt = open("%s/%s" % (self.stepdir, newname ), 'w')
            for head, seq in FastaParser("%s/%s" % (self.stepdir, f)):
                head = head.strip().split()[0]
                otpt.write("%s\t%s\n" % (head, id))
            otpt.close()
        if len(files)==0:
            self.message("No files to generate GROUPS...")          

class   MakeQualFile(DefaultStep):
    def __init__(self, PREV, q):
        ARGS =  {"qual": q}
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("FILE_qfile")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
        
    def performStep(self):
        f = self.find("fasta")[0]
        q = self.getInputValue("qual")
        self.message("Creating 'qual' file; '{0}' for sequences in {1}".format(q, f))
        newname = f.strip().split(".")[:-1]
        newname = "%s.qual" % (".".join(newname)) 
        otpt = open("%s/%s" % (self.stepdir, newname ), 'w')
        for head, seq in FastaParser("%s/%s" % (self.stepdir, f)):
            otpt.write(">%s\n" % (head))
            for k in seq:
                otpt.write("%s " % (q))
            otpt.write("\n")    
        otpt.close()            
        
class   AlignmentSummary(DefaultStep):
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("AlignmentSummary")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        self.project = __projectid__
        self.mailaddress = __email__
            
        f = self.find("fasta")[0] 
        
        ref =  self.getInputValue("ref")
        if ref == None:
            ref="e_coli2"
            
        th =  self.getInputValue("thresh")
        if th == None:
            th="0.1"    
        
        outfile_trim = "{}.altrimcoords".format(f)
        self.message("summarizing an alignment in %s" % (f) )
        k = "%spython %s/alignmentSummary.py -P %s -M %s -t 500 -p %s -i %s -o %s.alsum -T %s -x %s --output-trim %s" % (binpath, scriptspath, self.project, self.mailaddress, ref, f,f, th, binpath,outfile_trim)
        if YAPGlobals.debug_grid_tasks:
            k += " --debug-grid-tasks"
        if YAPGlobals.dummy_grid_tasks:
            k += " --dummy-grid-tasks"
        self.message(k)
        task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=True)
        task.wait()     
        
        x = loadLines(os.path.join(self.stepdir,outfile_trim))[-1].strip().split("\t")
        self.message("Potential trimming coordinates: %s - %s [peak = %s] [thresh = %s]" % (x[1], x[3], x[5], x[7]) )
        self.setOutputValue("trimstart", x[1])
        self.setOutputValue("trimend", x[3])
        self.setOutputValue("trimthresh", x[7])
        
        #self.failed = True
        

class   AlignmentPlot(DefaultStep):
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("AlignmentPlot")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        f = self.find("alsum")[0]
        
        ref =  self.getInputValue("ref")
        if ref == None:
            ref="e_coli"
            
        trimstart = self.getInputValue("trimstart")
        if trimstart==None:
            trimstart=0
        elif trimstart=="find":
            trimstart = self.find("trimstart")[0]
            
        trimend = self.getInputValue("trimend")
        if trimend == None:
            trimend=0
        elif trimend == "find":
            trimend = self.find("trimend")[0]
            
        trimthresh = self.getInputValue("trimthresh")
        if trimthresh == None:
            trimthresh=0
        elif trimthresh == "find":
            trimthresh = self.find("trimthresh")[0] 
            
        
        self.message("Adding trimmig marks at: %s - %s" % (trimstart, trimend))
        tmp = open("%s/alsum.r" % (self.stepdir), "w")
        tmp.write("source(\"%s/alignmentSummary.R\")\n" % (scriptspath))    
        tmp.write("batch2(\"%s\", ref=\"%s\", trimstart=%s, trimend=%s, thresh=%s )\n" % (f, ref, trimstart, trimend, trimthresh))
        tmp.close()
        k = "%sR CMD BATCH alsum.r" % (binpath)
        task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
        task.wait()
        
class   GroupRetriever(DefaultStep):
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("GroupRetriever")
        self.start()
        
    def performStep(self):
        minimum = self.getInputValue("mingroupmembers")
        if minimum==None:
            minimum=0
        
        group = self.find("group")[0]
        groups = defaultdict(int)
        otpt = open("{0}/{1}.groupstats".format(self.stepdir, group), "w")
        for line in loadLines("%s/%s" % (self.stepdir, group)):
            x = line.strip().split("\t")[1]
            groups[x]+=1
            
        keys  = sorted(groups, key=groups.get)
        keys.reverse()
        
        passinggroups=list()
        failinggroups = list()
        
        for k in keys:
            v = groups[k]
            if v>=minimum:
                flag="ok"
                passinggroups.append(k)
            else:
                flag="x"    
                failinggroups.append(k)
                
            self.message("{0:<25}:{1:>10}:{2}".format( k, v, flag))
            otpt.write("{0}\t{1}\t{2}\n".format(k,v, flag))
            
        otpt.close()
        
        if len(passinggroups)==0:
            self.message("There are not enough reads to analyze. See documentation for -g [currently set to {0}] and -x arguments.".format(minimum))
            self.failed=True
        
        if self.getInputValue("report") in [None, "passing"]:
            groupnames  = "-".join(passinggroups)   
        else:
            groupnames  = "-".join(failinggroups)   
        self.setOutputValue("groups", groupnames)

class   CDHIT_Preclust(DefaultStep):
    def __init__(self, nodeCPUs, ARGS, PREV):
        DefaultStep.__init__(self)
        if ARGS.has_key("T"):
            self.nodeCPUs = ARGS["T"]
        else:
            self.nodeCPUs=nodeCPUs
            ARGS["T"]=self.nodeCPUs

        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("CDHIT_Preclust")
        self.start()
                        
    def performStep(self):

        cdhitpath = globals()["cdhitpath"] ## to be found by locals()

        args = ""   

        arguments_copy = self.arguments.copy()
        n_parts = arguments_copy.pop("splits",10)
                
        for arg, val in arguments_copy.items():
            args = "%s -%s %s" % (args, arg, val) 

        fs = self.find("fasta")
        if len(fs)==0:
            fs.extend(self.find("mate1"))
            fs.extend(self.find("mate2"))

        scratch_dir = os.path.join(self.stepdir,"scratch")
        rmrf(scratch_dir)
        if not os.path.isdir(scratch_dir):
            os.makedirs(scratch_dir)
        scratch_dir = "scratch"

        for f in fs:
            f_base=os.path.join(scratch_dir,os.path.basename(f))
            work_base="{f_base}.db-uniq-split".format(**locals())
            db_uniq_clstr="{f_base}.db-uniq.clstr".format(**locals())
            
            k = "{cdhitpath}cd-hit-dup -i {f} -o {f_base}.db-uniq -m false -e 4".format(**locals())
            
            self.message(k)
            task_dup = GridTask(template=defaulttemplate, name=self.stepname+"_dup", command=k, cpu=1,  mem_per_cpu=31, dependson=list(), cwd = self.stepdir,
                    flag_completion = True)
            self.waitOnGridTasks([task_dup],fail_policy="ExitOnFirst")

            k = "{cdhitpath}cd-hit-div.pl {f_base}.db-uniq {work_base} {n_parts}".format(**locals())
            
            self.message(k)
            task_split = GridTask(template=defaulttemplate, name=self.stepname+"_div", command=k, cpu=1,  dependson=[], cwd = self.stepdir, flag_completion = True)
            self.waitOnGridTasks([task_split],fail_policy="ExitOnFirst")
            tasks_454 = []
            splits_fasta=[]
            splits_clstr=[]
            for i_div in range(n_parts):
                db_uniq_split="{}-{}".format(work_base,i_div)
                n_cpu_454=4
                
                k = """{cdhitpath}cd-hit-454 -i {db_uniq_split} -o {db_uniq_split}.nr1 -c 0.99 -n 10 -b 10 -g 0 -M 0 -aS 0.0 -T {n_cpu_454} && \
                       {cdhitpath}clstr_rev.pl {db_uniq_clstr} {db_uniq_split}.nr1.clstr > {db_uniq_split}.nr1-all.clstr && \
                       {cdhitpath}clstr_sort_by.pl < {db_uniq_split}.nr1-all.clstr > {db_uniq_split}.nr1-all.sort.clstr""".format(**locals())
                
                self.message(k)
                
                tasks_454.append(GridTask(template=defaulttemplate, name=self.stepname+"_div_454", 
                    command=k, cpu=n_cpu_454,  dependson=[], cwd = self.stepdir, flag_completion = True))
                
                splits_fasta.append("{db_uniq_split}.nr1".format(**locals()))
                splits_clstr.append("{db_uniq_split}.nr1-all.sort.clstr".format(**locals()))

            self.waitOnGridTasks(tasks_454,fail_policy="ExitOnFirst")
            
            ## concatenate results from splits
            
            k = """cat {} > {}.nr1-all.sort.clstr && \
                   cat {} > {}.nr1""".format(" ".join(splits_clstr),
                           work_base,
                           " ".join(splits_fasta),
                           work_base)
            
            self.message(k)
            
            task_cat = GridTask(template=defaulttemplate, name=self.stepname+"_cat", command=k, cpu=1,  dependson=[], cwd = self.stepdir, flag_completion = True)
            self.waitOnGridTasks([task_cat],fail_policy="ExitOnFirst")
            
            n_cpu_454=self.nodeCPUs
            
            k = """{cdhitpath}cd-hit-454 -i {work_base}.nr1 -o {work_base}.nr1.nr2 -c 0.98 -n 10 -b 10 -g 0 -M 0 -aS 0.0 -T {n_cpu_454} && \
                   {cdhitpath}clstr_rev.pl {work_base}.nr1-all.sort.clstr {work_base}.nr1.nr2.clstr > {work_base}.nr1.nr2-all.clstr && \
                   {cdhitpath}clstr_sort_by.pl < {work_base}.nr1.nr2-all.clstr > {work_base}.nr1.nr2-all.sort.clstr && \
                   mv {work_base}.nr1.nr2 {f}.cdhit && \
                   mv {work_base}.nr1.nr2-all.sort.clstr {f}.cdhit.clstr""".format(**locals())
            
            self.message(k)
            task_all_454 = GridTask(template=defaulttemplate, name=self.stepname+"_all_454", 
                    command=k, cpu=n_cpu_454,  dependson=[], cwd = self.stepdir, flag_completion = True)
            self.waitOnGridTasks([task_all_454],fail_policy="ExitOnFirst")
            
            if not YAPGlobals.debug_grid_tasks:
                GridTask(template=defaulttemplate, name=self.stepname+"_cleanup_454", 
                                            command="rm -rf scratch", cpu=1,  dependson=[], cwd = self.stepdir).wait()
            

class   CDHIT_454(DefaultStep):
    def __init__(self, nodeCPUs, ARGS, PREV):
        DefaultStep.__init__(self)
        if ARGS.has_key("T"):
            self.nodeCPUs = ARGS["T"]
        else:
            self.nodeCPUs=nodeCPUs
            ARGS["T"]=self.nodeCPUs     

        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("CDHIT_454")
        self.start()
                        
    def performStep(self):

        args = ""   
                
        for arg, val in self.arguments.items():
            args = "%s -%s %s" % (args, arg, val) 

        fs = self.find("fasta")
        if len(fs)==0:
            fs.extend(self.find("mate1"))
            fs.extend(self.find("mate2"))
        tasks=list()
        for f in fs:
            k ="%scd-hit-454 -i %s -o %s.cdhit %s" % (cdhitpath, f, f, args)
            self.message(k)
            task = GridTask(template=defaulttemplate, name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir)
#           if self.nodeCPUs>2:
#               task = GridTask(template=defaulttemplate, name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir, debug=True)
#           else:
#               task = GridTask(template="himem.q", name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir, debug=True)
            tasks.append(task)
        
        for task in tasks:
            task.wait() 
            
class   CDHIT_EST(DefaultStep):
    def __init__(self, nodeCPUs, ARGS, PREV):
        DefaultStep.__init__(self)
        
        if ARGS.has_key("T"):
            self.nodeCPUs = ARGS["T"]
        else:
            self.nodeCPUs=nodeCPUs
            ARGS["T"]=self.nodeCPUs     
        
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("CDHIT_EST")
        self.start()
                
    def performStep(self):
        f = self.find("fasta")[0]
        args = ""
        dist = 1
        for arg, val in self.arguments.items():
            args = "%s -%s %s" % (args, arg, val) 
            if arg == "c":
                dist = dist - (float(val)) 
        
        k ="%scd-hit-est -i  %s -o %s._%s_.cdhit %s" % (cdhitpath, f, f, dist, args)
        
        self.message(k)
        task = GridTask(template="pick", name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir)
#       if self.nodeCPUs>2:
#           task = GridTask(template="defaulttemplate", name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir, debug=True)
#       else:
#           task = GridTask(template="himem.q", name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir, debug=True)
        task.wait() 

class   CDHIT_Perls(DefaultStep):
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("CDHITperls")
        #self.nodeCPUs=nodeCPUs
        self.start()
        
    def performStep(self):
        x = self.find("cdhitclstr")
        tasks=list()
        for cluster in x:
            k = "%sclstr2tree.pl %s > %s.tre" % (cdhitpath, cluster, cluster)
            self.message(k)
            task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
            tasks.append(task)
            
            k = "%sclstr_size_histogram.pl %s > %s.hist.tab.txt " % (cdhitpath, cluster, cluster)
            self.message(k)
            task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
            tasks.append(task)
            
            k = "%sclstr_size_stat.pl %s  > %s.size.tab.txt" % (cdhitpath, cluster, cluster)
            self.message(k)
            task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
            tasks.append(task)
            
            
        for task in tasks:
            task.wait()                 
    
class   MakeAccnosFromName(DefaultStep):
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("MakeAccnosFromName")
        self.start()
        
    def performStep(self):
    
        n = self.find("name")  
        assert len(n) == 1
        n = n[0]
    
        k = "{}python {}MakeAccnosFromName.py -n {}  --min-cluster-size {}".format(binpath, 
                scriptspath, n, self.getInputValue("min_cluster_size"))    
        self.message(k)
        task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
        task.wait() 

class   CDHIT_Mothurize(DefaultStep):
    def __init__(self, ARGS, PREV):
        DefaultStep.__init__(self)
        #self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("CDHIT_Mothurise")
        self.start()
        
    def performStep(self):
    
        ### mode can be "name" or None (to produce a name file)
        ### mode can be "#.## (to produce a labeled ra/sa/list combo)
        m = self.getInputValue("mode")
        if m == None:
            m = "name"
        
        modeswitch = "-o %s" % (m)
        
        ### is there an optional names file?
        n = self.find("name")   
        if len(n)>0:
            n = n[0]
            nameswitch = "-n %s" % (n)
        else:
            nameswitch = "" 
    
        ### is there a required cluster file
        clst = self.find("cdhitclstr")
        
        if len(clst)>0:
            k = "%spython %sCDHIT_mothurize_clstr.py -c %s  %s %s" % (binpath, scriptspath, clst[0], nameswitch, modeswitch)    
            self.message(k)
            task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
            task.wait() 
            
        else:
            self.failed=True

class   R_defaultplots(DefaultStep):
    def __init__(self, INS, ARGS, PREV):
        DefaultStep.__init__(self)  
        self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("R_plots")
        self.start()
    def performStep(self):
        f = self.find("taxsummary")
        anno = self.find("annotation")[0]
        
        tasks = list()
    
        script = open("%s/script.r" % (self.stepdir), "w")
        script.write("""source("%sConsTaxonomyPlots.R")\n""" % (scriptspath))
        
        for file in f:
            dist = ".%s"% (self.getInputValue("dist"))
            
            if file.find(dist)>-1 and file.find("seq")>-1 :
                script.write("""makeDefaultBatchOfPlots("%s", "%s", fileprefix="SEQnum")\n""" % (anno, file))
    
                
            elif file.find(dist)>-1 and file.find("otu")>-1 :
                script.write("""makeDefaultBatchOfPlots("%s", "%s", fileprefix="OTUnum")\n""" % (anno, file))
            
        script.close()
        
        k = "%sR CMD BATCH script.r" % (binpath)    
        
        self.message(k)
        task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
        task.wait()
                
class   R_OTUplots(DefaultStep):
    def __init__(self, INS, ARGS, PREV):
        DefaultStep.__init__(self)  
        self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("R_plots_otu")
        self.start()
    def performStep(self):
    
        ####OTUS
        f = self.find("fasta")
        tasks = list()

        #script = open("%s/script.r" % (self.stepdir), "w")
        #script.write("""source("%sOtuReadPlots.r")\n""" % (scriptspath))
        
        for file in f:
            if file.find("annotated.fasta")>0:
                k = """grep ">" %s | awk '{FS = "|"; OFS="\t"} {print $4, $5}' > %s.otustats""" % (file, file)
                task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
                tasks.append(task)
                #script.write("""makeBatch("%s.otustats")\n""" % (file))
                
        ####COVERAGE
        f = self.find("clcassemblystats")
        
        #for file in f:
                #script.write("""makeBatchCoverage("%s")\n""" % (file))     
        
        #script.close()
        
        ### make sure all conversions are complete
        for task in tasks:
            task.wait()
        
        k = "%sR CMD BATCH %sOtuReadPlots.r" % (binpath, scriptspath)
        self.message(k)
        
        task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
        task.wait()             
                
class   R_rarefactions(DefaultStep):
    def __init__(self, INS, ARGS, PREV):
        DefaultStep.__init__(self)  
        self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("R_rarefactions")
        self.start()
    def performStep(self):
        for k in "r_nseqs,rarefaction,r_simpson,r_invsimpson,r_chao,r_shannon,r_shannoneven,r_coverage".strip().split(","):
            f = self.find(k)
        k = "%sR CMD BATCH %srarefactions.R" % (binpath, scriptspath)   
        self.message(k)
        task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
        task.wait()                 

class   AlignmentTrim(DefaultStep):
    def __init__(self, INS, ARGS, PREV):
        DefaultStep.__init__(self)
        self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("AlignmentTrim")
        self.start()
                        
    def performStep(self):
        f = self.find("fasta")[0]
        args = ""                   
        for arg, val in self.arguments.items():
            if val.startswith("find"):
                val=self.find(val.split(":")[1])[0]
            args = "%s -%s %s" % (args, arg, val) 
        
        k ="%spython %salignmentTrimmer.py %s -I %s" % (binpath, scriptspath, args, f)
        self.message(k)
        task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  dependson=list(), cwd = self.stepdir)
        task.wait() 
        
class   AnnotateClusters(DefaultStep):
    def __init__(self, INS, ARGS, PREV):
        DefaultStep.__init__(self)
        self.setInputs(INS)
        self.setArguments(ARGS)
        self.setPrevious(PREV)
        self.setName("Annotate")
        self.start()
                        
    def performStep(self):
        l = self.find("list")
        t = self.find("taxonomy")
        f = self.find("fasta")
        g = self.find("group")
        
        self.message(l)
        self.message(t)
        self.message(f)
        self.message(g)
            
        if len(l)==0 or len(t)==0 or len(f)==0 or len(g) == 0:
            self.failed=True
        else:
        
            tasks=list()    
            for fasta in f:
                dist = fasta.split("_")[-2]
                
                for tax in t:
                    if tax.find(dist)>-1 and tax.find("otu")==-1:
                        k = "%spython %sRetrieve.py %s %s %s %s %s" % (binpath, scriptspath, dist, l[0], tax, g[0], fasta)
                        self.message(k)
                        task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  dependson=list(), cwd = self.stepdir)
                        tasks.append(task)
            
            for task in tasks:
                task.wait() 

            
#################################################
##      Functions
##
    ################################################
    ### Read in a file and return a list of lines
    ###
def loadLines(x):
    try:
        fp = open(x, "r")
        cont=fp.readlines()
        fp.close()
        #print "%s line(s) loaded."  % (len(cont))
    except:
        cont=""
        #print "%s cannot be opened, does it exist? " % ( x )   
    return cont

def unlist(struct):
    for x in struct:
        if type(x) is tuple or type(x) is list or type(x) is set :
            for y in unlist(x):
                yield y
        else:
            yield x
    
def init(id, e, maxnodes = 250, update=0.1):
    global __projectid__
    global __email__
    global __admin__
    global BOH
    global MOTHUR
    global QS
    
    __projectid__ = id
    __email__ = e
    __admin__ = 'rsanka@jcvi.org'
    
    BOH = BufferedOutputHandler()
    MOTHUR = MothurCommandInfo(path=mothurpath)
    QS = TaskQueueStatus(update = update, maxnodes=maxnodes)
    
    return dict(BOH=BOH,QS=QS,MOTHUR=MOTHUR)
    
def revComp(string):
    global transtab
    string=string.upper()
    #reverse
    string = string [::-1]
    return string.translate(transtab)
    
#################################################
##      Arguments
##

#################################################
##      Begin
##

from string import maketrans
inttab=  "ACGTN"
outtab = "TGCAN"
transtab = maketrans(inttab, outtab)

pool_open_files = BoundedSemaphore(value=YAPGlobals.open_files_max, verbose=False)

rcfilepath = os.environ["YAP_RC"]

binpath = os.environ["YAP_DEPS"]+"/"
scriptspath = os.environ["YAP_SCRIPTS"]+"/"

mothurpath  = os.path.join(binpath,"mothur-current/")
cdhitpath   = os.path.join(binpath,"cdhit-current/")

print """\
        rcfilepath = {}\
        binpath = {}\
        scriptspath = {}\
        mothurpath = {}
        """.format(rcfilepath,binpath,scriptspath,mothurpath)

defaulttemplate = "default.q"


#################################################
##      Finish
#################################################

