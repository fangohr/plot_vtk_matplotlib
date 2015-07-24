import plot_vtk_matplotlib_2D as vtk_2d

# Load our skyrmion from a nanodisk
# of 80 nm diameter in the range [-40, 40]
testing_p = vtk_2d.plot_vtk_matplotlib('skyrmion.vtu', z_max=0.41, z_min=0.39)
testing_p.extract_data()

# Plot the z component
testing_p.plot_vtk(-40, 40, -40, 40, arrow_colour='k', save_fig='sk_mz.pdf')

# plot the x component and arrows in black
testing_p.plot_vtk(-40, 40, -40, 40, m_component='mx',
                   arrow_colour='k', save_fig='sk_mx_black.pdf')

# Quiver plot
testing_p.plot_quiver()
