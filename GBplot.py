import scipy.io as io
import plotly.offline as py
import plotly.graph_objs as go
import numpy as np
import os
import plotly.io as pio
py.init_notebook_mode(connected=True)

class Plot:

#########################################################################################
#  This class is only meant to plot all kinds of plots that the paper, the Coupling     #
#  Molecular Dynamics and Machine Learning for Predicting Aluminum Segregation to       #
#  Magnesium Grain Boundaries Required. It requires plotly, numpy and scipy libraries.  #
#  to be installed on your computer. Also 'data_Mg_GBperatom_seg_2Al_dump.mat', the     #
#  data file, should be at the same location with 'GBplot.py' doc.                      #
#########################################################################################
    
    # io = __import__('scipy.io')
    # py = __import__('plotly.offline')
    # go = __import__('plotly.graph_objs')
    # np = __import__('numpy')
    # os = __import__('os')
    
    def __init__(self):
        self.mat = io.loadmat('data_Mg_GBperatom_seg_2Al_dump.mat')
        self.length_A = self.mat['A'].shape[1]


    def peratom_energy_dist(self, gb_num, export = False):

    #########################################################################################
    #  Input:                                                                               #
    #        gb_num: The id of the Grain Boundary, ranging from 0 to 30.                    #
    #        export: Export .svg image if True, default is False.                          #
    #  Output:                                                                              #
    #        The plot of atom at their z axis position vs thier segregation energy.         #
    #########################################################################################
        
        i = gb_num 
        mat = self.mat

        #extract data from the specified Grain Boundary. 
        atom_pos = mat['A']['peratom'][0,i][0,0]['pos']
        
        segE = mat['A']['Eseg'][0,i]
        #check whether this is a valid data?
        n1 = segE[:,0] != 0 
        segE = segE[n1,:]
        
        atom_ID = segE[:,0].astype(int) - 1
        atom_zaxis = atom_pos[atom_ID,2]
        
        gbE = mat['A']['GBE'][0,i][0][0]
        gb_tilt = mat['A']['GB_tilt'][0,i][0]
        gb_norm1 = mat['A']['GB_norm1'][0,i][0]
        
        #draw scatter plot of the data
        main = go.Scatter(x = atom_zaxis-min(atom_zaxis)-20, 
                    y = segE[:,1], 
                    marker={'color':'red', 'symbol':'circle','size':3}, 
                    mode='markers',
                    name='Test'
                    )
        
        data = [main]
        
        layout = go.Layout(autosize = False, 
                    height = 800,
                    width = 800, 
                    xaxis={'title':'Distance from GB', 'zeroline':False},
                    yaxis = {'title':'Segeregation Energy(eV)','zeroline':False},
                    shapes =[{'type':'line', 
                                'x0':-25,
                                'x1':25,
                                'y0':np.mean(segE[:,1]),
                                'y1':np.mean(segE[:,1]),
                                'line':{'width':0.7
                                    
                                }
                                }],
                    annotations = [
                        dict(
                        x = 0.08,
                        y = 0.925,
                        xref = 'paper',
                        yref = 'paper',
                        text = '$\\text{GB Energy: %4.1f}\ mJ/m^2$'%(gbE),
                        showarrow = False,
                        font = dict(size = 10)
                        ),
                        dict(
                        x = 0.08,
                        y = 0.90,
                        xref = 'paper',
                        yref = 'paper',
                        text = '$\\text{GB Tilt: [%d %d %d %d]}$'%(gb_tilt[0],
                                                                gb_tilt[1],
                                                                gb_tilt[2],
                                                                gb_tilt[3]),
                        showarrow = False,
                        font = dict(size = 10)
                        ),
                        dict(
                        x = 0.08,
                        y = 0.875,
                        xref = 'paper',
                        yref = 'paper',
                        text = '$\\text{GB Normal: [%d %d %d %d]}$'%(gb_norm1[0],
                                                                    gb_norm1[1],
                                                                    gb_norm1[2],
                                                                    gb_norm1[3]
                        ),
                        showarrow = False,
                        font = dict(size = 10)
                        ),
                        dict(
                        x = 0.08,
                        y = 0.08,
                        xref = 'paper',
                        yref = 'paper',
                        text = '$\\text{Min Segregation: %4.2f eV}$'%(min(segE[:,1])),
                        showarrow = False,
                        font = dict(size = 10)
                        ),
                    ]
                    )
        fig = go.Figure(data = data, layout = layout)
        plot = py.iplot(fig)
        if export == True:
            if not os.path.exists('images'):
                os.mkdir('images')
            pio.write_image(fig, 'images/segE_dist%d.svg'%(gb_num))

        
    def pos_segE(self, gb_num, export = False):

    #########################################################################################
    #  Input:                                                                               #
    #        gb_num: The id of the Grain Boundary, ranging from 0 to 30.                    #
    #        export: Export .svg image if True, default is False.                          #
    #  Output:                                                                              #
    #        The plot of atom in their x, z position and color with thier energy level.The  #
    #        atom size has been scaled for better display.                                  #
    #########################################################################################

        i = gb_num
        mat = self.mat
        atom_pos = mat['A']['peratom'][0,i][0,0]['pos']
        num_atom = mat['A']['natoms'][0,i][0][0]
        segE = mat['A']['Eseg'][0,i]
        #check whether this is a valid data?
        n1 = segE[:,1] != 0 
        segE = segE[n1,:]
        atom_ID = segE[:,0].astype(int) - 1
        atom_zaxis = atom_pos[atom_ID,2]
        
        atomsize = 100*np.ptp(atom_pos[atom_ID,0])*np.ptp(atom_pos[atom_ID,2])/num_atom
        
        main = go.Scatter(x = atom_pos[atom_ID,0], 
                y = atom_pos[atom_ID,2],
                text = np.around(segE[:,1],3).tolist(),
                hoverinfo = 'text',
                marker={'size':atomsize/1.5, 'color':segE[:,1], 'colorscale':'Rainbow' ,
                'symbol':'circle','showscale':True}, 
                mode='markers'
                )
        data = [main]

        layout = go.Layout(autosize = False,
                    height = 800,
                    width = 1200/np.ptp(atom_pos[atom_ID,2])*np.ptp(atom_pos[atom_ID,0]),     
                    xaxis={'zeroline':False},
                    yaxis = {'zeroline':False}   
                    )
        fig = go.Figure(data = data, layout = layout)
        plot = py.iplot(fig)

        if export == True:
            if not os.path.exists('images'):
                os.mkdir('images')
            pio.write_image(fig, 'images/pos_segE%d.svg'%(gb_num))


    def correlation_segE_descriptor(self, descriptor_Name, subnum, export = False):
    ########################################################################################
    #  Input:                                                                              #     
    #        descriptor_name: The name of the physical descriptor. Some of the descriptor  #
    #        may not only have one term. The keywords are ['pe', 'pos', 'cna', 'stress',   #
    #        'centro_fnn', 'centro_snn', 'voronoi', 'f', 'coord']                          #
    #        subnum: subnum is the column variable of the descriptor we wish to plot.      #
    #        export: Export .svg image if True, default is False.                          # 
    #  Output:                                                                             #
    #        The plot of the segeregation energy value vs the descriptor value.            #
    ########################################################################################
        mat = self.mat
        for i in range(self.length_A):
            if i == 0:
                segE = mat['A']['Eseg'][0,0]
                #check whether this is a valid data?
                n1 = segE[:,0] != 0 
                segE = segE[n1,:]

                atom_ID = segE[:,0].astype(int) - 1
                descriptor_temp = mat['A']['peratom'][0,0][0,0][descriptor_Name][:,0]
                descriptor_all = descriptor_temp[atom_ID]
                segE_all = segE
            else:
                segE = mat['A']['Eseg'][0,i]
                #check whether this is a valid data?
                n1 = segE[:,0] != 0 
                segE = np.squeeze(segE[n1,:])

                atom_ID = segE[:,0].astype(int) - 1
                descriptor_temp = mat['A']['peratom'][0,i][0,0][descriptor_Name][:,0]
                descriptor_temp = descriptor_temp[atom_ID]
                descriptor_all = np.concatenate([descriptor_all, descriptor_temp])
                segE_all = np.concatenate([segE_all, segE])
        
        #calculate correlation coefficient
        R = np.corrcoef(np.squeeze(descriptor_all), segE_all[:,1])[0,1]
        
        #draw correlation map
        main = go.Scatter(x = np.squeeze(descriptor_all), 
                y = segE_all[:,1], 
                marker={'color':'red', 'symbol':'circle','size':3}, 
                mode='markers'
                )

        data = [main]

        layout = go.Layout(autosize = False, 
                height = 800,
                width = 800, 
                xaxis={'title':descriptor_Name, 'zeroline':False},
                yaxis = {'title':'Segregation Energy(eV)', 'zeroline':False},
                annotations = [
                        dict(
                        x = 0,
                        y = 1,
                        xref = 'paper',
                        yref = 'paper',
                        text = '$\\rho= \\text{%4.2f}$'%(R),
                        showarrow = False,
                        font = dict(size = 12)
                        )
                    ]
                )
        fig = go.Figure(data = data, layout = layout)
        plot = py.iplot(fig)

        if export == True:
            if not os.path.exists('images'):
                os.mkdir('images')
            pio.write_image(fig, 'images/Corr_segE_%s%d.svg'%(descriptor_Name, subnum))

    
    def correlation_map(self, absvalue = False, export = False):
    
    ########################################################################################
    #  Input:                                                                              #     
    #        absvalue: If True compute the absolute value of correlation coefficient map.  #
    #        export: Export .svg image if True, default is False.                          #
    #  Output:                                                                             #
    #        The plot of all correlation coefficients of all the feature including segE.   #
    ########################################################################################

        descriptor_keys = ['$pe$','$pos_{x}$','$pos_{y}$','$pos_{z}$','$CNA$','$Centro_{FNN}$',
                    '$Centro_{SNN}$','$Coord$','$f_{x}$','$f_{y}$','$f_{z}$','$\sigma_{x}$',
                    '$\sigma_{y}$','$\sigma_{z}$','$\sigma_{xy}$','$\sigma_{xz}$','$\sigma_{yz}$',
                    '$Vor_{vol}$','$Vor_{faces}$','$\sigma_{H}$','$f_{mag}$','$E_{seg}$']
        
        mat = self.mat
        
        for i in range(30):
            segE = mat['A']['Eseg'][0,i]
            #check whether this is a valid data?
            n1 = segE[:,0] != 0 
            segE = np.squeeze(segE[n1,:])
            atom_ID = segE[:,0].astype(int) - 1

            descriptor = mat['A']['peratom'][0,i][0,0]
            descriptor_temp = np.concatenate([descriptor['pos'],descriptor['pe'],descriptor['cna'],descriptor['centro_fnn'],
                                        descriptor['centro_snn'],descriptor['coord'],descriptor['f'],descriptor['stress'],
                                        descriptor['voronoi']], axis = 1)
            if i == 0:
                descriptor_all = descriptor_temp[atom_ID]
                segE_all = segE
            else:
                descriptor_temp = descriptor_temp[atom_ID]
                descriptor_all = np.concatenate([descriptor_all, descriptor_temp], axis = 0)
                segE_all = np.concatenate([segE_all, segE])

        descriptor_all[:,2] = abs(descriptor_all[:,2]-min(descriptor_all[:,2])-20)
        sigma_H = np.sum(descriptor_all[:,11:13], axis = 1)/3
        f_mag = np.linalg.norm(descriptor_all[:,8:11], axis = 1, ord = 2)

        feature = np.concatenate([descriptor_all, sigma_H[:,np.newaxis], f_mag[:,np.newaxis], segE_all[:,1][:,np.newaxis]], axis = 1)
        #calculate correlation coefficient
        R = np.corrcoef(feature, rowvar = False)
        
        #Optimal Leaf Algorithm for improving the order of the correlations
        import scipy.spatial.distance
        import scipy.cluster.hierarchy
        D = scipy.spatial.distance.pdist(R)
        tree = scipy.cluster.hierarchy.linkage(D, 'average')
        leafOrder = scipy.cluster.hierarchy.dendrogram(tree, no_plot = True)
        leafOrder = leafOrder['leaves']
        R = R[leafOrder,:]
        R = R[:,leafOrder]
        descriptor_keys = np.array(descriptor_keys)[leafOrder]
        if absvalue == True:
            R = abs(R)
            filename = 'images/AbsCorrelationMaps.svg'
        else:
            filename = 'images/CorrelationMaps.svg'
        trace = go.Heatmap(z=R,
                    x=descriptor_keys,
                    y=descriptor_keys,
                    colorscale = 'Portland')
        layout = dict(height = 1000, width = 1000)
        fig = go.Figure(data = [trace], layout = layout)
        if export == True:
            if not os.path.exists('images'):
                os.mkdir('images')
            pio.write_image(fig, filename)
        py.iplot(fig)