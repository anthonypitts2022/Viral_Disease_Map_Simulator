#author: Anthony Pitts

import random
import math
from matplotlib import pyplot as plt
import numpy as np

def normpdf(x, mean, sd):
    """
    Return the value of the normal distribution 
    with the specified mean and standard deviation (sd) at
    position x.
    """
    var = float(sd)**2
    denom = (2*math.pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

#returns probability of mortality
def pdeath(x, mean, sd):
    start = x-0.5
    end = x+0.5
    step =0.01    
    integral = 0.0
    while start<=end:
        integral += step * (normpdf(start,mean,sd) + normpdf(start+step,mean,sd)) / 2
        start += step            
    return integral    
    
recovery_time = 4 # recovery time in time-steps
virality = .6    # probability that a neighbor cell is infected in 
                  # each time step                                                  

class Cell(object):

    def __init__(self,x, y):
        self.x = x
        self.y = y 
        self.state = "S" # can be "S" (susceptible), "R" (resistant = dead), or 
                         # "I" (infected)
        self.time=0 #records the time
        
    def infect(self):
        self.state = "I" #"I" means "infected"
        self.time = 0 #records the time stamp of being infected
        
    def process(self, adjacent_cells):
        #this method is called to determine that change in state/time for a cell
        self.time= self.time + 1 
        if self.time>=recovery_time: #recovers cells that have been infected long enough
            self.state="S"
        if self.state=="I":
            if random.random()<=pdeath(self.time,.2,.2): #if cell should die
                self.state="R"
        if self.time>0 and self.state=="I":
            for ind_adj_cell in adjacent_cells:
                if ind_adj_cell.state=="S":
                    if (random.random()<=virality): #chance of infecting adj. cell
                        ind_adj_cell.infect()
        
class Map(object):
    
    def __init__(self):
        self.height = 150
        self.width = 150           
        self.cells = {} #key is tuple of (x,y) coordinates and values are Cell instances

    def add_cell(self, cell):
        #adds a cell to the map
        self.cells[(cell.x, cell.y)] = cell
        
    def display(self):
        image = np.zeros((self.width, self.height, 3))
        for i in range(0, self.width):
            for j in range(0, self.height):
                if ((i,j) in self.cells.keys()):
                    if self.cells[(i,j)].state=="S": #S makes cell green
                        image[i,j,0] = 0.0
                        image[i,j,1] = 1.0
                        image[i,j,2] = 0.0
                    elif self.cells[(i,j)].state=="R": #R makes cell green
                        image[i,j,0] = 0.0
                        image[i,j,1] = 1.0
                        image[i,j,2] = 0.0
                    elif self.cells[(i,j)].state=="I":#I makes cell red
                        image[i,j,0] = 1.0
                        image[i,j,1] = 0.0
                        image[i,j,2] = 0.0
                else: #coordinate not on map -- just make cell black
                    image[i,j,0] = 0.0
                    image[i,j,1] = 0.0
                    image[i,j,2] = 0.0
        plt.imshow(image)
    
    def adjacent_cells(self, x,y):
        adjacent_cells = [] #makes list of adjacent cells to each cell
        #right adjacent
        if ((x+1)<150):
            if (x+1,y) in self.cells.keys():
                adjacent_cells.append(self.cells[(x+1,y)])
        if ((x-1)>-1):
            if (x-1,y) in self.cells.keys():
                adjacent_cells.append(self.cells[(x-1,y)])
        if ((y+1)<150):
            if (x,y+1) in self.cells.keys():
                adjacent_cells.append(self.cells[(x,y+1)])
        if ((y-1)>-1):
            if (x,y-1) in self.cells.keys():
                adjacent_cells.append(self.cells[(x,y-1)])
        return adjacent_cells
    
    def time_step(self): #processes the increase in time
        for coordinate, eachCell in self.cells.items():
            eachCell.process(self.adjacent_cells(coordinate[0], coordinate[1]))
        self.display()
            
def read_map(filename):
    m = Map()
    with open(filename, 'r') as file:
        for line in file:
            #reads file for x y coordinates line by line
            line = line.replace("\n","")
            ind_cell_coordinates = line.split(",")
            m.add_cell(Cell(int(ind_cell_coordinates[0]), int(ind_cell_coordinates[1])))
    return m
