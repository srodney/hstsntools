# S.Rodney 2014.04.21
from bs4 import BeautifulSoup
import requests
import numpy as np

# BS4 doc:
# http://www.crummy.com/software/BeautifulSoup/bs4/doc

def checkProgram( pid, dayspan=10 ):
    """Fetch the visit status page, parse the visit info, print a report.
    """
    soup = fetchPID( pid )
    visdict = parseStatus(soup)
    print("--------------------------------------")
    print("HST Visit Status Update for program %i \n"%pid)
    mkReport(visdict, dayspan=dayspan )
    print("--------------------------------------")


def fetchPID( pid ):
    """ Read in the HST visit status page for the given program ID.
    """
    r  = requests.get("http://www.stsci.edu/cgi-bin/get-visit-status?id=%i&markupFormat=html&observatory=HST"%pid)
    data = r.text
    soup = BeautifulSoup(data)
    return( soup )

def parseStatus( soup ) :
    """ Parse the status page soup into a dictionary of visits.
    """
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

def mkReport( visdict, dayspan=8 ):
    """Construct a report of the visits that have been
    archived/executed/scheduled/withdrawn in
    the last/next day/week.
    """

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

    comingSoon = [ k for k in visdict.keys()
                          if ( visdict[k]['status'] in
                               ['Scheduled', 'Scheduling','Implementation']
                               and (0 < visdict[k]['daystostart'] < dayspan) )
                              or  (0 < visdict[k]['daystoend'] < dayspan)
                         ]

    print("Archived in the last %i days : "%dayspan +
          ','.join(archivedThisWeek).strip(','))
    print("Executed in the last %i days : "%dayspan +
          ','.join(executedThisWeek).strip(','))
    print("Scheduled for the last %i days : "%dayspan +
          ','.join(scheduledThisWeek).strip(','))
    print("Available in the next %i days : "%dayspan +
          ','.join(comingSoon).strip(','))


    datekey = lambda x : visdict[x]['enddate'].isoformat()

    doneVisits = sorted(doneLately,key=datekey,reverse=False)
    comingVisits = sorted(comingSoon,key=datekey,reverse=False)

    print('\n\n Visits Completed in the Last %i days:'%dayspan)
    for vis in doneVisits :
        datestr = visdict[vis]['enddate'].date().isoformat()
        timestr = visdict[vis]['enddate'].time().isoformat()[:5]
        target = visdict[vis]['targets']
        print( '%s %s %s  (%s)'%(vis,datestr,timestr,target))


    print('\n\n Visits Available for Scheduling in the Next %i days:'%dayspan)
    for vis in comingVisits :
        datestr = visdict[vis]['enddate'].date().isoformat()
        timestr = visdict[vis]['enddate'].time().isoformat()[:5]
        target = visdict[vis]['targets']
        print( '%s %s %s  (%s)'%(vis,datestr,timestr,target))



    # list archived/executed/withdrawn visits in the last day/week

    # list scheduled/withdrawn visits in the next day/week
    # parse plan windows
