//using System;
//using System.Collections.Generic;

//using Grasshopper.Kernel;
//using Rhino.Geometry;

//namespace DendroGH.Components
//{
//    public class MyComponent2_Sphere : GH_Component
//    {
//        /// <summary>
//        /// Initializes a new instance of the MyComponent1 class.
//        /// </summary>
//        public MyComponent2_Sphere()
//          : base("Create a Sphere From one point", "OneSphere",
//              "Create a sphere with center x,y,z and radius",
//              "Dendro", "Convert")
//        {
//        }

//        /// <summary>
//        /// Registers all the input parameters for this component.
//        /// </summary>
//        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
//        {
//            pManager.AddNumberParameter("radius_mm", "radius", "sphere radius", GH_ParamAccess.item);
//            pManager.AddNumberParameter("x coordinate of center", "center_x", "x coordinate of center", GH_ParamAccess.item);
//            pManager.AddNumberParameter("y coordinate of center", "center_y", "y coordinate of center", GH_ParamAccess.item);
//            pManager.AddNumberParameter("z coordinate of center", "center_z", "z coordinate of center", GH_ParamAccess.item);
           
//        }





//        /// <summary>
//        /// Registers all the output parameters for this component.
//        /// </summary>
//        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
//        {
//            pManager.AddGenericParameter("Sphere", "Sphere", "Create a sphere", GH_ParamAccess.item);
//        }

//        /// <summary>
//        /// This is the method that actually does the work.
//        /// </summary>
//        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
//        protected override void SolveInstance(IGH_DataAccess DA)
//        {
//            double radius = double.NaN;
//            double x = double.NaN;
//            double y = double.NaN;
//            double z = double.NaN;


//            if (!DA.GetData(0, ref radius)) return;
//            if (!DA.GetData(1, ref x)) return;
//            if (!DA.GetData(2, ref y)) return;
//            if (!DA.GetData(3, ref z)) return;



//            DendroVolume sphere = new DendroVolume(x,y,z,radius);
           

//            DA.SetData(0, new VolumeGOO(sphere));
//        }

//        /// <summary>
//        /// Provides an Icon for the component.
//        /// </summary>
//        protected override System.Drawing.Bitmap Icon
//        {
//            get
//            {
//                //You can add image files to your project resources and access them like this:
//                // return Resources.IconForThisComponent;
//                return null;
//            }
//        }

//        /// <summary>
//        /// Gets the unique ID for this component. Do not change this ID after release.
//        /// </summary>
//        public override Guid ComponentGuid
//        {
//            get { return new Guid("0f8e63aa-8265-4fa4-86e2-9565f105f767"); }
//        }
//    }
//}
