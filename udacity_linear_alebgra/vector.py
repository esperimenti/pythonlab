from math import sqrt,acos,pi
import numpy as np
from decimal import Decimal,getcontext

getcontext().prec = 30

class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            #self.coordinates = tuple(coordinates)
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates
    
    def plus(self,v):
        result = [a+b for a,b in zip(self.coordinates,v.coordinates)]
        return Vector(result)

    def minus(self,v):
        result = [a-b for a,b in zip(self.coordinates,v.coordinates)]
        return Vector(result)

    
    def times_scalar(self, c):
        result = [c*x for x in self.coordinates]
        return Vector(result)
    
    def magnitude(self):
        #v1
        #power_res_list = [x**2 for x in self.coordinates]
        #z = 0
        #for y in power_res_list:
        #    z = z + y
        #result = math.sqrt(z)
        #v2
        power_res_list = [x**2 for x in self.coordinates]
        result = sqrt(sum(power_res_list))
        return Decimal(result)
    
    def normalized(self):
        #return self.times_scalar(1/self.magnitude());
        return self.times_scalar(Decimal('1.0')/self.magnitude());
    
    
    def dot(self,v):
        '''Calculate inner product of 2 vectors'''
        #v1
        #return sum(np.array(self.coordinates)*np.array(v.coordinates))
        #v2
        return sum([x*y for x,y in zip(self.coordinates,v.coordinates)])
    
    def angle_with(self, vector):
        """
        Calculate the angle between two Vectors in Radians
        Where Theta == arccos(v1 dot v2 / |v1| dot |v2|)
        """
        vector_dot_prod = self.dot(vector)
        if self.magnitude() == 0 or vector.magnitude() == 0:
            raise ZeroDivisionError('A zero vector has no angle')

        magnitude_dot_prod = self.magnitude() * vector.magnitude()
        cos_angle = min(1, max(vector_dot_prod / magnitude_dot_prod, -1))
        return acos(cos_angle)
    
    #def angle_with(self,v):
        #v1
        #in_prod =  self.inner_product(v)
        #product_mag = self.magnitude()*v.magnitude()
        #return acos(in_prod/product_mag)
        #v2
    #    u1 = self.normalized()
    #    u2 = v.normalized()
    #    print('u1:',u1)
    #    print('u2:',u2)
    #    return acos(u1.dot(u2))
        
    def angle_degree(self,v):
        return self.angle_with(v)*180/pi
        
    
    def is_parallel_to(self,v):
        return (self.is_zero() or
                v.is_zero() or
                self.angle_with(v) in [0,pi])
        
    def is_zero(self,tolerance=1e-10):
        return self.magnitude() < tolerance
    
    def is_orthogonal_to(self,v,tolerance=1e-10):
        if abs(self.dot(v))< tolerance:        
            return True
        else:
            return False
    
    #def project_to(self, v):
    #    """
    #    Project one Vector onto another
    #    Where the resultant projected_vector = (v1 dot v2_norm) * v2norm
    #    """
    #    v_norm = v.normalized()
    #    return v_norm.times_scalar(self.dot(v_norm))
    
    def component_orthogonal_to(self,basis):
        try:
            projection = self.component_parallel_to(basis)
            return self.minus(projection)
        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e
                
    def component_parallel_to(self,basis):
        try:
            u = basis.normalized()
            weight = self.dot(u)
            return u.times_scalar(weight)
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e
                
    def cross_product(self,v):
        x1 = self.coordinates[0]
        y1 = self.coordinates[1]
        z1 = self.coordinates[2]
        x2 = v.coordinates[0]
        y2 = v.coordinates[1]
        z2 = v.coordinates[2]
        
        p1 = y1*z2-y2*z1
        p2 = -(x1*z2-x2*z1)
        p3 = x1*y2-x2*y1
        
        return Vector([p1,p2,p3])
     
    def area_parallelogram(self,v):
        #return sqrt(sum([x**2 for x in self.cross_product(v).coordinates]))
        return self.cross_product(v).magnitude()
    
    def area_triangle(self,v):
        return self.area_parallelogram(v)/2
