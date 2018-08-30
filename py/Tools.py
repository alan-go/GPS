import math
def llh_xyz(lat,long,h):
    long = long*math.pi/180
    lat  = lat*math.pi/180
    f = 1/298.257223563
    e2 = f*(2-f)
    a = 6378137.0
    N = a/math.sqrt(1-e2*math.sin(lat)*math.sin(lat))
    x = (N+h)*math.cos(lat)*math.cos(long)
    y = (N+h)*math.cos(lat)*math.sin(long)
    z = (N*(1-e2)+h)*math.sin(lat)
    return x,y,z


if __name__=="__main__":
    lat = 39.97928617
    long = 116.33892268
    h = 54.922
    x,y,z = llh_xyz(lat,long,h)
    print(x,y,z)
