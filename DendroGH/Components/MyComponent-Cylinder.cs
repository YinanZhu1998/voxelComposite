using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace DendroGH.Components
{
    public class MyComponent_Cylinder : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent_Cylinder class.
        /// </summary>
        public MyComponent_Cylinder()
          : base("Cylinder", "Cylinder",
              "Create a cylinder by radius and two points",
              "Dendro", "Convert")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("radius_mm", "radius", "sphere radius", GH_ParamAccess.item);
            pManager.AddNumberParameter("x coordinate of center1", "center1_x", "x coordinate of center", GH_ParamAccess.item);
            pManager.AddNumberParameter("y coordinate of center1", "center1_y", "y coordinate of center", GH_ParamAccess.item);
            pManager.AddNumberParameter("z coordinate of center1", "center1_z", "z coordinate of center", GH_ParamAccess.item);
            pManager.AddNumberParameter("x coordinate of center2", "center2_x", "x coordinate of center", GH_ParamAccess.item);
            pManager.AddNumberParameter("y coordinate of center2", "center2_y", "y coordinate of center", GH_ParamAccess.item);
            pManager.AddNumberParameter("z coordinate of center2", "center2_z", "z coordinate of center", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Cylinder", "Cylinder", "Create a cylinder", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double radius = 0;
            double x1 = double.NaN;
            double y1 = double.NaN;
            double z1 = double.NaN;
            double x2 = double.NaN;
            double y2 = double.NaN;
            double z2= double.NaN;


            if (!DA.GetData(0, ref radius)) return;
            if (!DA.GetData(1, ref x1)) return;
            if (!DA.GetData(2, ref y1)) return;
            if (!DA.GetData(3, ref z1)) return;
            if (!DA.GetData(4, ref x2)) return;
            if (!DA.GetData(5, ref y2)) return;
            if (!DA.GetData(6, ref z2)) return;



            DendroVolume cylinder = new DendroVolume();
            cylinder.CreateCylinderFromTwoPoint(radius,x1,y1,z1,x2,y2,z2);


            DA.SetData(0, new VolumeGOO(cylinder));
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
            get { return new Guid("78b0ff26-74cc-4057-b9fa-f336bcdf0c40"); }
        }
    }
}