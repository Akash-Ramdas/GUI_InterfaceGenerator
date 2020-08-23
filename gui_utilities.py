#Widget Libraries
from IPython.display import display
import ipywidgets as widgets
from ipywidgets import Button, Layout, Text, IntSlider, Output, VBox, HBox, Style, GridBox, Label
#Basic Libraries
import numpy as np

#pymatgen libraries
import pymatgen
from pymatgen import MPRester
from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#a = MPRester("OClHROzM9CGykIgeA")

#ase libraries
from ase import Atoms
from ase.build import surface
from ase.build import cut
from ase.build import bulk
from ase import cell
from ase.io import write
from ase.build import stack

#3D Viewer (nglview) libraries
import nglview
#Function to create all widgets used in GUI

def init_gui(a):

	def create_widgets():
		global title, material_1, material_2, view_3d, output_1, h1_slide, k1_slide, l1_slide;
		global h2_slide, k2_slide, l2_slide, list_interfaces, output_2, SID1_slide, SID2_slide;
		global gen_slabs, substrate_slide, gap_slide, merge_slabs, rep_x1, rep_y1, rep_z1, rep_x2;
		global rep_y2, rep_z2, rep_x3, rep_y3, rep_z3, rep_interface, rep_structures, down_interface;
		global down_structures, file_format;
		#Create all widgets here
		#Text as Buttons
		title = Button(
		    description="Interface Generator",
		    layout=Layout(height='auto', width='auto', grid_area= "title"),
		)

		#Input Texts
		material_1 = Text(
		    value='mp-149',
		    placeholder='Enter MPID (e.g: mp-2858) or Chemical Formula(e.g:ZrO2)',
		    description='Material 1:',
		    disabled=False,
		    layout=Layout(
		        height='auto', width='auto',
		        grid_area="mat_1",
		        #justify_items = "center"
		    )
		)

		material_2 = Text(
		    value='mp-352',
		    placeholder='Enter MPID (e.g: mp-2858) or Chemical Formula(e.g:ZrO2)',
		    description='Material 2:',
		    disabled=False,
		    layout=Layout(
		        height='auto', width='auto', 
		        grid_area="mat_2"
		        #justify_items = "center"
		    )
		)

		#Buttons
		view_3d = Button(
		    description="View 3D structure",
		    layout=Layout(height='auto', width='auto', grid_area="view_3d"),
		    tooltip='Click to view 3D structures'
		)

		list_interfaces = Button(
		    description="List Possible Slabs", 
		    disabled=True,
		    layout=Layout(height='auto', width='auto', grid_area="list_interfaces"),
		    tooltip='Click to list the possible slabs'
		)

		gen_slabs = Button(
		    description="Show Slabs", 
		    disabled=True,
		    layout=Layout(height='auto', width='auto', grid_area="gen_slabs"),
		    tooltip='Click to generate the desired slabs'
		)

		merge_slabs = Button(
		    description="Generate Interface", 
		    disabled=True,
		    layout=Layout(height='auto', width='auto', grid_area="merge_slabs"),
		    tooltip='Click to generate the desired interface'
		)

		rep_structures = Button(
		    description="Replicate Structures", 
		    disabled=True,
		    layout=Layout(height='auto', width='auto', grid_area="rep_structures"),
		    tooltip='Replicate the two structures'
		)

		rep_interface = Button(
		    description="Replicate Interface", 
		    disabled=True,
		    layout=Layout(height='auto', width='auto', grid_area="rep_interface"),
		    tooltip='Replicate the generated interface'
		)


		down_structures = Button(
		    description="Download Structures", 
		    disabled=True,
		    layout=Layout(height='auto', width='auto', grid_area="down_structures"),
		    tooltip='Download the two structures in desired formate'
		)

		down_interface = Button(
		    description="Download Interface", 
		    disabled=True,
		    layout=Layout(height='auto', width='auto', grid_area="down_interface"),
		    tooltip='Download the generated interface in desired format'
		)


		#Sliders
		h1_slide = widgets.IntSlider(
		    value=1, 
		    description='h_1', 
		    max=3, 
		    min=0,
		    layout=Layout(height='auto', width='auto', grid_area="h1_slide")
		);

		k1_slide = widgets.IntSlider(
		    value=1, 
		    description='k_1', 
		    max=3, 
		    min=0,
		    layout=Layout(height='auto', width='auto', grid_area="k1_slide")
		);
		l1_slide = widgets.IntSlider(
		    value=1, 
		    description='l_1', 
		    max=3, 
		    min=0,
		    layout=Layout(height='auto', width='auto', grid_area="l1_slide")
		);
		h2_slide = widgets.IntSlider(
		    value=1, 
		    description='h_2', 
		    max=3, 
		    min=0,
		    layout=Layout(height='auto', width='auto', grid_area="h2_slide")
		);

		k2_slide = widgets.IntSlider(
		    value=1, 
		    description='k_2', 
		    max=3, 
		    min=0,
		    layout=Layout(height='auto', width='auto', grid_area="k2_slide")
		);
		l2_slide = widgets.IntSlider(
		    value=1, 
		    description='l_2', 
		    max=3, 
		    min=0,
		    layout=Layout(height='auto', width='auto', grid_area="l2_slide")
		);

		SID1_slide = widgets.IntSlider(
		    value=0, 
		    description='SID_1', 
		    max=0, 
		    min=0,
		    layout=Layout(height='auto', width='auto', grid_area="SID1_slide")
		);
		SID2_slide = widgets.IntSlider(
		    value=0, 
		    description='SID_2', 
		    max=0, 
		    min=0,
		    layout=Layout(height='auto', width='auto', grid_area="SID2_slide")
		);

		substrate_slide = widgets.IntSlider(
		    value=1, 
		    description='Substrate (1 for material 1, 2 for material 2)', 
		    max=2, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="substrate_slide")
		);
		gap_slide = widgets.FloatSlider(
		    value=0.1, 
		    description= 'Gap between 2 interfaces desired in Angstrom', 
		    min = 0.1, 
		    max = 5, 
		    step = 0.1,
		    layout=Layout(height='auto', width='auto', grid_area="gap_slide")
		);

		#Int Text
		rep_x1 = widgets.BoundedIntText(
		    value=1, 
		    description='x1', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_x1")
		);

		rep_y1 = widgets.BoundedIntText(
		    value=1, 
		    description='y1', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_y1")
		);

		rep_z1 = widgets.BoundedIntText(
		    value=1, 
		    description='z1', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_z1")
		);

		rep_x2 = widgets.BoundedIntText(
		    value=1, 
		    description='x2', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_x2")
		);

		rep_y2 = widgets.BoundedIntText(
		    value=1, 
		    description='y2', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_y2")
		);

		rep_z2 = widgets.BoundedIntText(
		    value=1, 
		    description='z2', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_z2", )
		);


		rep_x3 = widgets.BoundedIntText(
		    value=1, 
		    description='x', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_x3")
		);

		rep_y3 = widgets.BoundedIntText(
		    value=1, 
		    description='y', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_y3")
		);

		rep_z3 = widgets.BoundedIntText(
		    value=1, 
		    description='z', 
		    max=10, 
		    min=1,
		    layout=Layout(height='auto', width='auto', grid_area="rep_z3")
		);
		#Dropdowns
		file_format = widgets.Dropdown(options=[".cif","POSCAR (VASP)", "Quantum Espresso input", ".xyz"],
		                               layout=Layout(height='auto', width='auto', grid_area="file_format"))
		#Outputs
		output_1 = widgets.Output(layout=Layout(height='auto', width='auto', grid_area= "output_1"));
		output_2 = widgets.Output(layout=Layout(height='auto', width='auto', grid_area= "output_2"));

	#Functionality of Widgets here
	#Utilities to convert pymatgen structures to ase structures for modification
	def build_ase_struc_3d(pymatgen_struc):
	    s = Atoms(pymatgen_struc.formula.replace(" ", ""),
	               positions = pymatgen_struc.cart_coords,
	               cell=pymatgen_struc.lattice.get_cartesian_coords(1),
	               pbc=[1,1,1]);
	    return s;
	def build_ase_struc_2d(pymatgen_struc):
	    s = Atoms(pymatgen_struc.formula.replace(" ", ""),
	               positions = pymatgen_struc.cart_coords,
	               cell=pymatgen_struc.lattice.get_cartesian_coords(1),
	               pbc=[1,1,0]);
	    return s;

	#Fetching materials data from MP database
	def get_3D_Structures(b):
	    global view, view_2, struc1,struc2, s1, s2, pressed;
	    pressed = 0;
	    struc1 = a.get_structure_by_material_id(material_1.value);
	    struc1 = SpacegroupAnalyzer(struc1).get_conventional_standard_structure()
	    s1 = build_ase_struc_3d(struc1)
	    view = nglview.show_ase(s1)  
	    view.add_representation('unitcell')
	    view.background = 'black'
	    view.camera = "orthographic"
	    
	    
	    
	    struc2 = a.get_structure_by_material_id(material_2.value);
	    struc2 = SpacegroupAnalyzer(struc2).get_conventional_standard_structure()
	    s2 = build_ase_struc_3d(struc2)
	    view_2 = nglview.show_ase(s2)  
	    view_2.add_representation('unitcell')
	    view_2.background = 'black'
	    view_2.camera = "orthographic"

	    with output_1:
	        output_1.clear_output()
	        display(view)
	        display(view_2)
	    
	    list_interfaces.disabled = False;
	    rep_structures.disabled = False;
	    down_structures.disabled = False;

	#Generating list possible slab terminations using slab generator from pymatgen    
	def list_possible_2D_Slabs(b):
	    global view, view_2, struc1,struc2, s1, s2, slabgen1, slabgen2, pressed;
	    pressed=0;
	    slabgen1 = SlabGenerator(struc1, (int(h1_slide.value),int(k1_slide.value),int(l1_slide.value)), 10, 10,center_slab=True);
	    slabgen2 = SlabGenerator(struc2, (int(h2_slide.value),int(k2_slide.value),int(l2_slide.value)), 10, 10,center_slab=True)
	    
	    with output_2:
	        output_2.clear_output()
	        print("For Material 1")
	        print("There are actually now %s terminations that can be generated" %(len(slabgen1.get_slabs())))
	        i = 0;
	        print("S-ID \t Miller Indices Polar? Symmetric? ")
	        for slab in slabgen1.get_slabs():
	            print(i, "\t", slab.miller_index,"\t", slab.is_polar(), "\t", slab.is_symmetric())
	            i = i+1;
	        i = 0;
	        SID1_slide.max = len(slabgen1.get_slabs())-1
	        
	        print("For Material 2")
	        print("There are actually now %s terminations that can be generated" %(len(slabgen2.get_slabs())))
	        print("S-ID \t Miller Indices  Polar? Symmetric? ")
	        for slab in slabgen2.get_slabs():
	            print(i, "\t", slab.miller_index,"\t", slab.is_polar(), "\t", slab.is_symmetric())
	            i = i+1;
	            
	        SID2_slide.max = len(slabgen2.get_slabs())-1
	        gen_slabs.disabled = False

	#ASE functionality to merge two slabs of the material
	def merge_2D_slabs(b):
	    global view, view_2, struc1,struc2, s1, s2, slabgen1, slabgen2,s3, view_3, pressed;
	    pressed = 1;
	    a1_tmp = (s2.cell.lengths()[0]/(s1.cell.lengths()[0]));
	    a2_tmp = (s2.cell.lengths()[1]/(s1.cell.lengths()[1]));
	    a1 = 1;
	    b1 = 1;
	    a2 = 1;
	    b2 = 1;
	    if(a1_tmp >= 1):
	        a1 = int(np.round(a1));
	    else:
	        b1 = int(np.round(1/a1));
	    
	    if(a2_tmp >= 1):
	        a2 = int(np.round(a1));
	    else:
	        b2 = int(np.round(1/a1));
	        
	    s1 = s1.repeat((a1,a2,1));
	    s2 = s2.repeat((b1,b2,1));
	    
	    if(substrate_slide.value == 1):
	        cell_final = s1.cell.copy();
	        fix_val = 0;
	    else:
	        cell_final = s2.cell.copy();
	        fix_val = 1;
	    
	    
	    vac = gap_slide.value/2;
	    s1.center(vacuum=vac, axis=2)
	    s2.center(vacuum=vac, axis=2)
	    
	    s3, s1_strain, s2_strain = stack(s1, s2, axis= 2, fix=0, cell = cell_final, maxstrain=None, distance=None, output_strained = True);
	    
	    strain1 = np.sqrt(((s1_strain.cell - s1.cell).sum(axis=0)**2).sum());
	    strain2 = np.sqrt(((s2_strain.cell - s2.cell).sum(axis=0)**2).sum());
	    
	    view_3 = nglview.show_ase(s3)  
	    view_3.add_representation('unitcell')
	    view_3.background = 'black'
	    view_3.camera = "orthographic"
	    
	    with output_1:
	        output_1.clear_output()
	        print("Strain on cell of Material 1 is", strain1);
	        print("Strain on cell of Material 2 is", strain2);
	        display(view_3)

	#Generate desired 2D slabs from the list of slabs generated by pymatgen
	def generate_2D_Slabs(b):
	    global view, view_2, struc1,struc2, s1, s2, slabgen1, slabgen2, pressed;
	    pressed = 0;
	    struc1 = slabgen1.get_slabs()[SID1_slide.value].get_sorted_structure();
	    s1 = build_ase_struc_2d(struc1);
	    
	    view = nglview.show_ase(s1)  
	    view.add_representation('unitcell')
	    view.background = 'black'
	    view.camera = "orthographic"

	    struc2 = slabgen2.get_slabs()[SID2_slide.value].get_sorted_structure();
	    s2 = build_ase_struc_2d(struc2);

	    
	    view_2 = nglview.show_ase(s2)  
	    view_2.add_representation('unitcell')
	    view_2.background = 'black'
	    view_2.camera = "orthographic"
	    
	    
	    with output_1:
	        output_1.clear_output()
	        display(view)
	        display(view_2)
	    
	    merge_slabs.disabled = False;
	    rep_interface.disabled = False;
	    down_interface.disabled = False;

	#Function to replicate slabs/3D models in x,y,z
	def replicate_structures(b):
	    global view, view_2, s1, s2;
	    s1 = s1.repeat((rep_x1.value,rep_y1.value,rep_z1.value));
	    s2 = s2.repeat((rep_x2.value,rep_y2.value,rep_z2.value));
	    view = nglview.show_ase(s1)  
	    view.add_representation('unitcell')
	    view.background = 'black'
	    view.camera = "orthographic"
	    view_2 = nglview.show_ase(s2)  
	    view_2.add_representation('unitcell')
	    view_2.background = 'black'
	    view_2.camera = "orthographic"
	    with output_1:
	        output_1.clear_output()
	        display(view)
	        display(view_2)

	#Function to replicate interface in x,y,z
	def replicate_interface(b):
	    global view_3, s3;
	    s3 = s3.repeat((rep_x3.value,rep_y3.value,rep_z3.value));
	    view_3 = nglview.show_ase(s3)  
	    view_3.add_representation('unitcell')
	    view_3.background = 'black'
	    view_3.camera = "orthographic"
	    with output_1:
	        output_1.clear_output()
	        display(view_3)

	#Function to download structures
	def download_structures(b):
	    global s1, s2;
	    if(file_format.value == "POSCAR (VASP)"):
	        write("struc1.POSCAR",s1, format="vasp" );
	        write("struc2.POSCAR",s2, format="vasp" );
	    elif(file_format.value == ".cif"):
	        write("struc1.cif",s1, format="cif" );
	        write("struc2.cif",s2, format="cif" );
	    elif(file_format.value == ".xyz"):
	        write("struc1.xyz",s1, format="xyz" );
	        write("struc2.xyz",s2, format="xyz" );
	    elif(file_format.value == "Quantum Espresso input"):
	        write("struc1.qe.in",s1, format="espresso-in" );
	        write("struc2.qe.in",s2, format="espresso-in" );

	#Function to download interface
	def download_interface(b):
	    global s3;
	    if(file_format.value == "POSCAR (VASP)"):
	        write("struc3.POSCAR",s3, format="vasp" );
	    elif(file_format.value == ".cif"):
	        write("struc3.cif",s1, format="cif" );
	    elif(file_format.value == ".xyz"):
	        write("struc3.xyz",s3, format="xyz" );
	    elif(file_format.value == "Quantum Espresso input"):
	        write("struc3.qe.in",s3, format="espresso-in" );

	create_widgets()
	view_3d.on_click(get_3D_Structures);
	list_interfaces.on_click(list_possible_2D_Slabs)
	gen_slabs.on_click(generate_2D_Slabs)
	merge_slabs.on_click(merge_2D_slabs)
	rep_structures.on_click(replicate_structures)
	rep_interface.on_click(replicate_interface)
	down_structures.on_click(download_structures)
	down_interface.on_click(download_interface)
	grid=GridBox(children=[title, material_1, material_2, view_3d, output_1, h1_slide, k1_slide, l1_slide,
	                      h2_slide, k2_slide, l2_slide, list_interfaces, output_2, SID1_slide, SID2_slide,
	                      gen_slabs, substrate_slide, gap_slide, merge_slabs, rep_x1, rep_y1, rep_z1, rep_x2,
	                      rep_y2, rep_z2, rep_x3, rep_y3, rep_z3, rep_interface, rep_structures, down_interface,
	                      down_structures, file_format],
	        layout=Layout(
	        width = '100%',
	        grid_template_columns='10% 12.5% 12.5% 10% 28% 13.5% 13.5%',
	        grid_template_rows='auto auto auto auto auto auto auto auto auto auto suto auto auto',
	        #grid_gap='0px',
	        grid_template_areas='''
	        "title title title title title title title"
	        ". mat_1 mat_1 . output_1 . ."
	        ". mat_2  mat_2 . output_1 rep_x1 rep_x2"
	        ". view_3d view_3d . output_1 rep_y1 rep_y2"
	        "h1_slide h1_slide h2_slide h2_slide output_1 rep_z1 rep_z2"
	        "k1_slide k1_slide k2_slide k2_slide output_1 rep_structures rep_structures"
	        "l1_slide l1_slide l2_slide l2_slide output_1 rep_x3 rep_x3"
	        ". list_interfaces list_interfaces . output_1 rep_y3 rep_y3"
	        ". . . . output_1 rep_z3 rep_z3"
	        ". . . . output_1 rep_interface rep_interface"
	        "output_2 output_2 output_2 output_2 output_1 . ."
	        "SID1_slide SID1_slide SID2_slide SID2_slide output_1 file_format file_format"
	        ". gen_slabs gen_slabs . output_1 down_structures down_structures"
	        "substrate_slide substrate_slide gap_slide gap_slide output_1 down_interface down_interface"
	        ". merge_slabs merge_slabs . output_1 . ."
	        ''',
	        #justify_items='center'
	        )
	       )
	display(grid)
