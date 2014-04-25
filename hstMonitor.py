#! /usr/bin/env python
# Author: S.Rodney
# Created : 2014.04.21

import sys
try:
    from apscheduler.scheduler import Scheduler
except ImportError :
    print("Error:  hstMonitor requires APScheduler.")
    print("        http://pythonhosted.org/APScheduler")
    print("Install it via pip (or, if you prefer, easy_install)")
    print("   pip install apscheduler")
    sys.exit()


sched = Scheduler( standalone=True )

@sched.cron_schedule( minute='*' )
def testReportMinutes( ):
    """Print the time every minute (for testing the cron-like scheduling and
    job persistence).
    """
    import time
    print( time.asctime() )

@sched.cron_schedule( day='*', hour=0, minute=0, second=0 )
def checkProgramDaily( pid, emailto='',emailuser='',emailpass='',
                       verbose=False ):
    """Every day at 00:00:00, fetch the visit status page,
    parse the visit info, print a report. If the email options
    are provided, the report is also emailed.
    """
    import time
    soup = fetchPID( pid )
    visdict = parseStatus(soup)
    preface = """
%s
Daily HST Visit Status Update for program %i

"""%( time.asctime(), pid )

    footer = """
Visit Status page: http://www.stsci.edu/cgi-bin/get-visit-status?id=%i

MAST archive page: http://archive.stsci.edu/hst/search.php?sci_pep_id=%i&action=Search&outputformat=HTML_Table&max_records=100
--------------------------------------
"""%( pid, pid )
    report = reportDone(visdict, 1 )
    if report and verbose :
        print( preface + report + footer )
    elif verbose :
        print( "Daily HST Visit status: nothing to report.")

    if report and emailto and emailuser and emailpass :
        # send a notice for visits completed in the last 1 day
        sendgmail( emailuser, emailpass, emailto,
                   'HST Visits completed : PID %i'%pid, preface+report )

@sched.cron_schedule( day_of_week='sun')
def checkProgramWeekly( pid, emailto='',emailuser='',emailpass='',
                        verbose=False ):

    """Each Sunday at 00:00:00, fetch the visit status page,
    parse the visit info, and print a report.  If the email options
    are provided, the report is also emailed.
    """
    import time
    dayspan=7
    soup = fetchPID( pid )
    visdict = parseStatus(soup)

    preface = """
%s
Weekly HST Visit Status Update for program %i

"""%( time.asctime(), pid )

    footer = """
    Visit Status page: http://www.stsci.edu/cgi-bin/get-visit-status?id=%i

MAST archive page: http://archive.stsci.edu/hst/search.php?sci_pep_id=%i&action=Search&outputformat=HTML_Table&max_records=100
--------------------------------------
"""%( pid , pid )

    report = reportComing(visdict, 7 )
    if report and verbose :
        print( preface + report + footer )
    elif verbose :
        print( "Daily HST Visit status: nothing to report.")

    if report and emailto and emailuser and emailpass :
        # send a notice for visits scheduled in the next 7 days
        sendgmail( emailuser, emailpass, emailto,
                   'Weekly HST schedule update for PID %i'%pid, preface + report + footer )


def sendgmail( username, password, toaddr, subject, message, ccaddr=''):
    """Send an email using gmail.
    """
    import smtplib
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEText import MIMEText

    fromaddr = "%s@gmail.com"%username
    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['CC'] = ccaddr
    msg['Subject'] = subject
    msg.attach( MIMEText( message, 'plain'))

    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login( username, password)
    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()


def fetchPID( pid ):
    """ Read in the HST visit status page for the given program ID.
    """
    try:
        from bs4 import BeautifulSoup
    except ImportError :
        print("Error:  hstMonitor requires BeautifulSoup4.")
        print("        http://www.crummy.com/software/BeautifulSoup")
        print("Install it via pip (or, if you prefer, easy_install)")
        print("   pip install beautifulsoup4")
        sys.exit()

    try:
        import requests
    except ImportError :
        print("Error:  hstMonitor requires .")
        print("        http://docs.python-requests.org/en/latest")
        print("Install it via pip  (or, if you prefer, easy_install)")
        print("   pip install requests")
        sys.exit()

    r  = requests.get("http://www.stsci.edu/cgi-bin/get-visit-status?id=%i&markupFormat=html&observatory=HST"%pid)
    data = r.text
    soup = BeautifulSoup(data)
    return( soup )

def parseStatus( soup ) :
    """ Parse the status page soup into a dictionary of visits.
    """
    import numpy as np
    from dateutil.parser import parse as dtparse
    # from dateutil.relativedelta import relativedelta
    import datetime
    today = datetime.datetime.today()

    visdict = {}

    trowlist = soup.find_all('tr')
    for trow in trowlist :
        rowtext = trow.getText().strip().split('\n')

        rowtext = np.array( rowtext )[ np.where( [ len(rt)  for rt in rowtext ] )[0] ]

        if rowtext[0] == 'Visit' : continue
        visit=rowtext[0]
        status=rowtext[1]
        targets=rowtext[2].replace('ANY','')
        configs=rowtext[3]

        startdate = dtparse( '9999-09-09' )
        enddate = dtparse( '9999-09-09')

        if status in ['Archived', 'Executed', 'Scheduled'] :
            startdate = dtparse( rowtext[4].split('-')[0] )
            enddate = dtparse( rowtext[5].split('-')[0] )
        elif status in ['Scheduling','Implementation'] :
            planwindowlist = rowtext[4:]
            for planwindow in planwindowlist :
                if '-' not in planwindow : continue
                try :
                    startdate0 = dtparse( planwindow.split('-')[0] )
                    enddate0 = dtparse( planwindow.split('-')[1].split('(')[0] )
                except : 
                    continue
                daystostart0 = (startdate0 - today ).days
                daystoend0 = (startdate0 - today ).days
                if daystostart0>=0 and daystoend0>=0 :
                    startdate = startdate0
                    enddate = enddate0
                    break
        daystostart =  ( startdate - today ).days
        daystoend =  ( enddate - today ).days

        visdict[visit] = {'visit':visit,
                          'status':status,
                          'targets':targets,
                          'configs':configs,
                          'startdate':startdate,
                          'enddate':enddate,
                          'daystostart':daystostart,
                          'daystoend':daystoend,
                          }
    return( visdict )


def reportDone( visdict, dayspan=8 ):
    """Construct a report of the visits that have been
    archived/executed/scheduled in the last <dayspan> days, and return
    it to the user.
    If nothing has been done lately, returns 0.
    """
    daynames = ['Mon','Tue','Wed','Thu','Fri','Sat','Sun']

    archivedThisWeek = [ k for k in visdict.keys()
                         if ( visdict[k]['status']=='Archived'
                              and -dayspan < visdict[k]['daystoend'] < 1 )
                         ]
    scheduledThisWeek = [ k for k in visdict.keys()
                          if ( visdict[k]['status']=='Scheduled'
                              and -dayspan < visdict[k]['daystostart'] < 1 )
                         ]
    executedThisWeek = [ k for k in visdict.keys()
                          if ( visdict[k]['status']=='Executed'
                              and -dayspan < visdict[k]['daystoend'] < 1 )
                         ]

    doneLately = [ k for k in visdict.keys()
                         if ( visdict[k]['status'] in
                              [ 'Archived', 'Executed', 'Scheduled' ]
                              and -dayspan < visdict[k]['daystoend'] < 1 )
                         ]

    if len(doneLately)==0 : return('')

    report = "\n"
    report += "Archived in the last %i days : "%dayspan
    report += ','.join(archivedThisWeek).strip(',') + '\n'
    report += "Executed in the last %i days : "%dayspan
    report += ','.join(executedThisWeek).strip(',') + '\n'
    report += "Scheduled for the last %i days : "%dayspan
    report += ','.join(scheduledThisWeek).strip(',') + '\n'

    datekey = lambda x : visdict[x]['enddate'].isoformat()
    doneVisits = sorted(doneLately,key=datekey,reverse=False)

    report += '\n\n Visits Completed in the Last %i days:\n'%dayspan
    for vis in doneVisits :
        datestr = visdict[vis]['enddate'].date().isoformat()
        timestr = visdict[vis]['enddate'].time().isoformat()[:5]
        target = visdict[vis]['targets']
        weekday = daynames[ visdict[vis]['enddate'].date().weekday() ]
        report += '%s   %s %s %s  (%s)\n'%(
            vis, weekday, datestr, timestr, target )

    report += '\n'
    return( report )


def reportComing( visdict, dayspan=8 ):
    """Construct a report of the visits that have been
    scheduled to be done in the next <dayspan> days.
    If nothing is scheduled, returns an empty string.
    """
    daynames = ['Mon','Tue','Wed','Thu','Fri','Sat','Sun']

    comingSoon = [ k for k in visdict.keys()
                          if ( visdict[k]['status'] in
                               ['Scheduled', 'Scheduling','Implementation']
                               and (0 < visdict[k]['daystostart'] < dayspan) )
                              or  (0 < visdict[k]['daystoend'] < dayspan)
                         ]

    datekey = lambda x : visdict[x]['enddate'].isoformat()
    comingVisits = sorted(comingSoon,key=datekey,reverse=False)

    if len(comingVisits)==0 : return('')

    report = '\n\n Visits Scheduled (or schedulable) for the Next %i days:\n'%dayspan
    for vis in comingVisits :
        datestr = visdict[vis]['enddate'].date().isoformat()
        timestr = visdict[vis]['enddate'].time().isoformat()[:5]
        target = visdict[vis]['targets']
        weekday = daynames[ visdict[vis]['enddate'].date().weekday() ]
        report += '%s   %s %s %s  (%s)\n'%(
            vis, weekday, datestr,timestr,target)

    report += '\n'
    return( report )


def mkReport( visdict, lookback=1, lookahead=7 ):
    """Construct a report of the visits that have been
    archived/executed in the last <lookback> days
    and those scheduled for execution in the next <lookahead> days.
    """
    report1 = reportDone( visdict, dayspan=lookback )
    report2 = reportComing( visdict, dayspan=lookahead )
    return( report1 + report2 )



if __name__ == "__main__":
    import argparse
    import datetime

    parser = argparse.ArgumentParser(
        description='Fetch the visit status page, parse the visit info, print a report.')

    # Required positional argument
    parser.add_argument('PID', type=str, help='HST Program ID to check.')
    pidlist = PID.split(',')


    # optional arguments
    # parser.add_argument('--lookback', metavar='N', type=float,
    #                     help='Number of days before today to search for completed visits.',
    #                     default=1)
    # parser.add_argument('--lookahead', metavar='N', type=float,
    #                     help='Number of days after today to search for scheduled visits.',
    #                     default=8)
    parser.add_argument('--quiet',  action='store_true', help='Suppress all stdout print statements', default='False')

    mailpar = parser.add_argument_group( "Options for e-mailing reports via gmail")
    parser.add_argument('--emailto', metavar='A,B,C', type=str, help='email addresses to send reports to.', default='')
    parser.add_argument('--emailuser', metavar='X', type=str, help='gmail username, for sending reports.', default='')
    parser.add_argument('--emailpass', metavar='Y', type=str, help='gmail password, for sending reports.', default='')

    argv = parser.parse_args()
    verbose = not argv.quiet

    for pid in pidlist :
        # every day, check for newly-completed visits
        checkProgramDaily( pid, emailto=argv.emailto,

                           emailuser=argv.emailuser, emailpass=argv.emailpass,
                           verbose=verbose )

        # every week, check for newly-scheduled visits
        checkProgramWeekly( pid, emailto=argv.emailto,
                            emailuser=argv.emailuser, emailpass=argv.emailpass,
                            verbose=verbose  )

    sched.start()



