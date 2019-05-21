using UnityEngine;

namespace SimpleKeplerOrbits
{

	[System.Serializable]
	public class AttractorData
	{
		public Transform AttractorObject;
        //massa do sol
		public float AttractorMass = 1000;
        //constante gravitacional
		public float GravityConstant = 0.1f;
	}
}