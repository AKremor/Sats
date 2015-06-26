#!/Python27/python
print 'Content-type: text/html'
print

__author__ = 'Anthony'
import satTLE
import cgi
import json

formData = cgi.FieldStorage()
if formData['calling'].value == 'getSatellites':
    satelliteList = satTLE.loadTLE()
    names = dict((key, value) for (key, value) in satelliteList)
    print json.dumps(names)