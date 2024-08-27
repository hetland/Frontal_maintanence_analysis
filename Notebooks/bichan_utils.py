import matplotlib.pyplot as plt
import cmocean.cm as cmo
import numpy as np
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
proj = ccrs.LambertConformal(central_longitude=-92, central_latitude=29)
coast_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='k', facecolor='0.8')


def to_rho(var, grid=None):
    if grid is None:
        grid = var.grid
    if var.dims[-1] != 'xi_rho':
        var = grid.interp(var, 'X', to='center', boundary='extend')
    if var.dims[-2] != 'eta_rho':
        var = grid.interp(var, 'Y', to='center', boundary='extend')
    return var


def add_derivatives(ds, grid, q='salt'):
    
    qs = ds[q]
    
    #################
    # Flow properties
    
    u = to_rho(ds.u, grid)
    v = to_rho(ds.v, grid)
    ds['ke'] = 0.5*(u**2 + v**2)
    
    #############################
    # Flow and property gradients
    
    ds['dqdx'] = to_rho(grid.derivative(qs, 'X'), grid)    # defined at rho-points
    ds['dqdy'] = to_rho(grid.derivative(qs, 'Y'), grid)    # defined at rho-points
    
    ds['dudx'] = grid.derivative(ds.u, 'X', boundary='extend')  # defined at rho-points
    ds['dvdy'] = grid.derivative(ds.v, 'Y', boundary='extend')  # defined at rho-points
    ds['dvdx'] = to_rho(grid.derivative(ds.v, 'X', boundary='extend'), grid)  # defined at rho-points
    ds['dudy'] = to_rho(grid.derivative(ds.u, 'Y', boundary='extend'), grid)  # defined at rho-points
    
    ###########################
    # Invariant flow properties
    
    # Vorticity:  v_x - u_y
    ds['zeta'] = (ds.dvdx - ds.dudy)/ds.f
    ds['zeta'].name = 'Normalized vorticity'

    # Divergence: u_x + v_y
    ds['delta'] = (ds.dudx + ds.dvdy)/ds.f
    ds['delta'].name = 'Normalized divergence'

    # Major axis of deformation
    ds['sigma'] = ( np.sqrt( (ds.dudx-ds.dvdy)**2 + (ds.dvdx+ds.dudy)**2 ) )/ds.f
    ds['sigma'].name = 'Normalized total strain'

    ##################################
    # Principle deformation components

    ds['lminor'] = 0.5 * (ds.delta - ds.sigma)
    ds['lminor'].name = 'lambda minor'

    ds['lmajor'] = 0.5 * (ds.delta + ds.sigma)
    ds['lmajor'].name = 'lambda major'
    
    #############################################
    # Along- and cross-frontal velocity gradients
    
    # angle is wrt x, so need to do arctan2(y, x)
    ds['phi_cf'] = np.arctan2(ds.dqdy, ds.dqdx)
    ds['phi_af'] = ds.phi_cf + np.pi/2.0

    ds['du_cf'] = ( ds.dudx*np.cos(ds.phi_cf)**2 + ds.dvdy*np.sin(ds.phi_cf)**2 
               + (ds.dudy + ds.dvdx)*np.sin(ds.phi_cf)*np.cos(ds.phi_cf) )/ds.f

    ds['du_af'] = ( ds.dudx*np.cos(ds.phi_af)**2 + ds.dvdy*np.sin(ds.phi_af)**2
              + (ds.dudy + ds.dvdx)*np.sin(ds.phi_af)*np.cos(ds.phi_af) )/ds.f
    
    ######################################################################
    # Strain efficiency, and the angle of the strain relative to the front
    
    ds['tilde_sigma_n'] = ds.du_cf - ds.du_af
    ds['tilde_sigma_n'].name = 'Normalized effective strain'
    
    ds['strain_efficiency'] = ds.tilde_sigma_n / ds.sigma
    ds['strain_efficiency'].name = 'Strain efficiency'
    
    ############################
    # The frontogenesis function
    
    # Dimensional frontogenesis function
    Dgradq_i = - ds.dudx*ds.dqdx - ds.dvdx*ds.dqdy
    Dgradq_j = - ds.dudy*ds.dqdx - ds.dvdy*ds.dqdy
    ds['Ddelq2'] = (ds.dqdx*Dgradq_i + ds.dqdy*Dgradq_j)
    ds['Ddelq2'].name = 'Frontogenesis function'

    # Density gradients squared
    ds['gradq2'] = ds.dqdx**2 + ds.dqdy**2
    ds['gradq2'].name = r'$(\nabla q)^2$'

    # Normalized frontogenesis function
    ds['nFGF'] = ds.Ddelq2 / (ds.gradq2 * ds.f)
    ds['nFGF'].name = r'Normalized Frontogenesis Function'
    
    ds['Area'] = 1.0/(ds.pn * ds.pm)
    
    return ds

def anim(ds, time_idxs, 
         anim_name='test.mp4'):

    if not os.path.exists('./frames/'):
        os.makedirs('./frames/')

    frame = 0
    
    for n in time_idxs:

        fig, axs = plt.subplots(1, 4, figsize=(12, 5))
        
        ax = axs[0];
        pc = ds.salt.isel(ocean_time=n).plot(x='x_rho', y='y_rho', ax=ax, cmap=cmo.haline,
                                             cbar_kwargs={'label': 'Salinity [g/kg]'})
        pc.set_clim(15, 55)
        ax.set_aspect(1.0)
        ax.set_xlabel('x [km]')
        ax.set_ylabel('y [km]')
        ax.set_xticklabels([str(x/1000.0) for x in ax.get_xticks()])
        ax.set_yticklabels([str(y/1000.0) for y in ax.get_yticks()])
        ax.set_title('Salinity')
        tstring = str(ds.ocean_time[n].values)
        day = tstring[8:10]
        hour = tstring[11:13]
        ax.text(0.05, 0.02, 
                '{hour} h {day} d\n$\\Delta x$ = 1 km'.format(day=day, hour=hour) ,
                horizontalalignment='left',
                verticalalignment='bottom', color='w',
                transform=ax.transAxes)
        
        ax = axs[1]
        pc = ds.zeta.isel(ocean_time=n).plot(x='x_rho', y='y_rho', ax=ax, cmap=cmo.curl,
                                             cbar_kwargs={'label': r'$\zeta/f$'})
        pc.set_clim(-5, 5)
        ax.set_aspect(1.0)
        ax.set_xlabel('x [km]')
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.set_xticklabels([str(x/1000.0) for x in ax.get_xticks()])
        ax.set_yticklabels([str(y/1000.0) for y in ax.get_yticks()])
        ax.set_title('Vorticity')
        
        ax = axs[2]
        pc = ds.delta.isel(ocean_time=n).plot(x='x_rho', y='y_rho', ax=ax, cmap=cmo.balance, 
                                             cbar_kwargs={'label': r'$\delta/f$'})
        pc.set_clim(-5, 5)
        ax.set_aspect(1.0)
        ax.set_xlabel('x [km]')
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.set_xticklabels([str(x/1000.0) for x in ax.get_xticks()])
        ax.set_yticklabels([str(y/1000.0) for y in ax.get_yticks()])
        ax.set_title('Divergence')
        
        ax = axs[3]
        gradq2 = (ds.dqdx**2 + ds.dqdy**2)
        gradq2.name = r'log_{10}($\nabla$ salt)$^2)$'
        pc = np.log10(gradq2).isel(ocean_time=n).plot(x='x_rho', y='y_rho', ax=ax, cmap=cmo.amp,
                                             cbar_kwargs={'label': r'log$_{10}$($\nabla$ salt)$^2)$'})
        pc.set_clim(-8, -5)
        ax.set_aspect(1.0)
        ax.set_xlabel('x [km]')
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.set_xticklabels([str(x/1000.0) for x in ax.get_xticks()])
        ax.set_yticklabels([str(y/1000.0) for y in ax.get_yticks()])
        ax.set_title(r'Salt gradient')

        fig.savefig('frames/frame_%04d.png' % frame, facecolor='w')
        frame += 1
        
        plt.close(fig)

    os.system('/global/homes/h/hetland/.conda/envs/seahorce/bin/ffmpeg -y -r 10 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 15 ' + anim_name)
    
    
def add_box(ds, ax, ihslice, draw_labels=False, timestr=None):
    ax.add_feature(coast_10m)

    ax.plot(ds.lon_rho[ihslice['eta_rho'].start, ihslice['xi_rho']],
            ds.lat_rho[ihslice['eta_rho'].start, ihslice['xi_rho']],
            '-', color='k', lw=2, transform=ccrs.PlateCarree())
    ax.plot(ds.lon_rho[ihslice['eta_rho'].stop, ihslice['xi_rho']],
            ds.lat_rho[ihslice['eta_rho'].stop, ihslice['xi_rho']],
            '-', color='k', lw=2, transform=ccrs.PlateCarree())
    ax.plot(ds.lon_rho[ihslice['eta_rho'], ihslice['xi_rho'].start],
            ds.lat_rho[ihslice['eta_rho'], ihslice['xi_rho'].start],
            '-', color='k', lw=2, transform=ccrs.PlateCarree())
    ax.plot(ds.lon_rho[ihslice['eta_rho'], ihslice['xi_rho'].stop],
            ds.lat_rho[ihslice['eta_rho'], ihslice['xi_rho'].stop],
            '-', color='k', lw=2, transform=ccrs.PlateCarree())
    
    ax.set_extent([-94.9, -91.7, 27.7, 29.9], crs=ccrs.PlateCarree())
    ax.set_title('')
    
    if timestr:
        ax.text(0.97, 0.03, timestr, ha='right', va='bottom', fontsize=16,
        transform=ax.transAxes)
    
    gl = ax.gridlines(draw_labels=draw_labels, x_inline=False, y_inline=False, 
                      xlocs=np.arange(-94.5,-90, 0.5), ylocs=np.arange(27, 31, 0.5))

    # manipulate `gridliner` object to change locations of labels
    gl.top_labels = False
    gl.right_labels = False