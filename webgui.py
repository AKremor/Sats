import web
from math import radians
import satTLE
from web import form
urls = ('/', 'index')
render = web.template.render('templates/')



sat_form = form.Form(
    #form.Textbox('satname', description = 'Satellite'),
    form.Dropdown('satname', sorted(list(satTLE.loadTLE().keys())), description='Satellite'),

    form.Button('submit', type='submit', description='Submit'))
    
class index:
    def GET(self):
        f = sat_form()
        satname = ''
        az = ''
        el = ''
        return render.index(f,satname,az,el)
        
    def POST(self):
        f = sat_form()
        satname =  web.input().satname
        obsvLLA = [radians(-36.377518), radians(145.400044), 100]
        az, el = satTLE.returnAzEl(satname, obsvLLA)
        return render.index(f,satname,az,el)
        
if __name__ == "__main__":
    app = web.application(urls,globals())
    app.run()
    
    
    
    
    
    
    

