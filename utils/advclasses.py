class Mesh():
    """
    Mesh class, for keeping all the mesh data.
    Input: directory [string] -> place where the domain files are.
    """
    def __init__(self, directory=""):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            domain_x = np.loadtxt(directory+"domain_x.dat")
        except IOError:
            print "IOError with domain_x.dat"
        try:
            #We have to avoid ghost cells!
            domain_y = np.loadtxt(directory+"domain_y.dat")[3:-3]
        except IOError:
            print "IOError with domain_y.dat"
        self.xm = domain_x #X-Edge
        self.ym = domain_y #Y-Edge
        
        self.xmed = 0.5*(domain_x[:-1] + domain_x[1:]) #X-Center
        self.ymed = 0.5*(domain_y[:-1] + domain_y[1:]) #Y-Center
        
        #(Surfaces taken from the edges)
        #First we make 2D arrays for x & y, that are (theta,r)
        T,R = meshgrid(self.xm, self.ym)
        R2  = R*R
        self.surf = 0.5*(T[:-1,1:]-T[:-1,:-1])*(R2[1:,:-1]-R2[:-1,:-1])

class Parameters():
    """
    Class for reading the simulation parameters.
    input: string -> name of the parfile, normally variables.par
    """
    def __init__(self, directory=''):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            params = open(directory+"variables.par",'r') #Opening the parfile
        except IOError:                  # Error checker.
            print  paramfile + " not found."
            return
        lines = params.readlines()     # Reading the parfile
        params.close()                 # Closing the parfile
        par = {}                       # Allocating a dictionary
        for line in lines:             #Iterating over the parfile
            name, value = line.split() #Spliting the name and the value (first blank)
            try:
                float(value)           # First trying with float
            except ValueError:         # If it is not float
                try:
                    int(value)         #                   we try with integer
                except ValueError:     # If it is not integer, we know it is string
                    value = '"' + value + '"'
            par[name] = value          # Filling the dictory
        self._params = par             # A control atribute, actually not used, good for debbuging
        for name in par:               # Iterating over the dictionary
            exec("self."+name.lower()+"="+par[name]) #Making the atributes at runtime


class Field(Mesh, Parameters):
    """
    Field class, it stores all the mesh, parameters and scalar data 
    for a scalar field.
    Input: field [string] -> filename of the field
           staggered='c' [string] -> staggered direction of the field. 
                                      Possible values: 'x', 'y', 'xy', 'yx'
           directory='' [string] -> where filename is
           dtype='float64' (numpy dtype) -> 'float64', 'float32', 
                                             depends if FARGO_OPT+=-DFLOAT is activated
    """
    def __init__(self, field, staggered='c', directory='', dtype='float64'):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        Mesh.__init__(self, directory) #All the Mesh attributes inside Field!
        Parameters.__init__(self, directory) #All the Parameters attributes inside Field!

        #Now, the staggering:
        if staggered.count('x')>0:
            self.x = self.xm[:-1] #Do not dump the last element
        else:
            self.x = self.xmed
        if staggered.count('y')>0:
            self.y = self.ym[:-1]
        else:
            self.y = self.ymed
    
        self.data = self.__open_field(directory+field,dtype) #The scalar data is here.

    def __open_field(self, f, dtype):
        """
        Reading the data
        """
        field = fromfile(f, dtype=dtype)
        return field.reshape(self.ny, self.nx)

    def plot(self, log=False, cartesian=False, cmap='Oranges_r', **karg):
        """
        A layer to plt.imshow or pcolormesh function.
        if cartesian = True, pcolormesh is launched.
        """
        ax = gca()
        if log:
            data = np.log(self.data)
        else:
            data = self.data
        if cartesian:
            T,R = meshgrid(self.x,self.y)
            X = R*cos(T)
            Y = R*sin(T)
            pcolormesh(X,Y,data,cmap=cmap,**karg)
        else:
            ax.imshow(data, cmap = cmap, origin='lower',aspect='auto',
                      extent=[self.x[0],self.x[-1],self.y[0],self.y[-1]],
                      **karg)
            
    def contour(self, log=False, cartesian=False, **karg):
        if log:
            data = np.log(self.data)
        else:
            data = self.data
        ax = gca()
        T,R = meshgrid(self.x,self.y)
        if cartesian:
            X = R*cos(T)
            Y = R*sin(T)
            contour(X,Y,data,**karg)
        else:
            contour(T,R,data,**karg)


