using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace DendroGH.Components
{
    public class afterGraduationFibers : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent_cylinders class.
        /// </summary>
        public afterGraduationFibers()
          : base("fibers", "fibers",
              "Create voxel based fibers",
              "Dendro", "Convert")
        {
            
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("aspectRatioMean", "aspectRatioMean", "aspectRatioMean", GH_ParamAccess.item);
            pManager.AddNumberParameter("aspectRatioDe", "aspectRatioDe", "aspectRatioDe", GH_ParamAccess.item);
            pManager.AddNumberParameter("diameterMean", "diameterMean", "diameterMean", GH_ParamAccess.item);
            pManager.AddNumberParameter("diameterDe", "diameterDe", "diameterDe", GH_ParamAccess.item);
            pManager.AddNumberParameter("orientationMean", "orientationMean", "orientationMean", GH_ParamAccess.item);
            pManager.AddNumberParameter("orientationDe", "orientationDe", "orientationDe", GH_ParamAccess.item);
            pManager.AddNumberParameter("volumeFraction", "volumeFraction", "volumeFraction", GH_ParamAccess.item);
        }
            

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("fibers", "fibers", "fibers", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double arMean = double.NaN;
            double arDe = double.NaN;
            double diameterMean = double.NaN;
            double diameterDe = double.NaN;
            double orientationMean = double.NaN;
            double orientationDe = double.NaN;
            double volumeFraction = 0;

            if (!DA.GetData(0, ref arMean)) return;
            if (!DA.GetData(1, ref arDe)) return;
            if (!DA.GetData(2, ref diameterMean)) return;
            if (!DA.GetData(3, ref diameterDe)) return;
            if (!DA.GetData(4, ref orientationMean)) return;
            if (!DA.GetData(5, ref orientationDe)) return;
            if (!DA.GetData(6, ref volumeFraction)) return;

            DendroVolume fibers = new DendroVolume();
            fibers.FibersGenerator(arMean, arDe, diameterMean, diameterDe, orientationMean, orientationDe, volumeFraction);

            DA.SetData(0, new VolumeGOO(fibers));
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
            get { return new Guid("319a5eb3-c54e-4093-8df9-86e4a13d0dbe"); }
        }
    }
}