using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace DendroGH.Components
{
    public class MyComponent_cylinders : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent_cylinders class.
        /// </summary>
        public MyComponent_cylinders()
          : base("Cylinders", "RanCylinders",
              "Create a cylinder by radius and two points",
              "Dendro", "Convert")
        {
            
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("length of bounding box", "length_mm", "length_mm", GH_ParamAccess.item);
            pManager.AddNumberParameter("width of bounding box", "width_mm", "width_mm", GH_ParamAccess.item);
            pManager.AddNumberParameter("height of bounding box", "height_mm", "height_mm", GH_ParamAccess.item);
           
            pManager.AddNumberParameter("number of points", "num", "num", GH_ParamAccess.item);
        }
            

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Cylinders", "Cylinders", "Create cylinders", GH_ParamAccess.item);
        }
        
        

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
        double length = double.NaN;
        double width = double.NaN;
        double height = double.NaN;
 
        double Num = 0;


        if (!DA.GetData(0, ref length)) return;
        if (!DA.GetData(1, ref width)) return;
        if (!DA.GetData(2, ref height)) return;
        
        if (!DA.GetData(3, ref Num)) return;




        DendroVolume cylinders = new DendroVolume();
        cylinders.CreateCylinderFromRandomPoint(length, width, height, Num);


        DA.SetData(0, new VolumeGOO(cylinders));
    }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("318a5eb3-c54e-4093-8df9-86e4a13d0dbe"); }
        }
    }
}