using UnityEngine;
using System.Collections;

public class RevealNightEarthTexture : MonoBehaviour {

	Transform tfLight;

	void Start () {
		var goLight = GameObject.Find ("RevealingLight");
		if(goLight)
		{
			tfLight = goLight.transform;
		}
	
	}
	
	void Update () {
	
		if(tfLight)
		{
            tfLight.position = new Vector3(transform.position.x, transform.position.y, transform.position.z + 4.3f);
			GetComponent<Renderer>().material.SetVector("_LightPos", tfLight.position);
			GetComponent<Renderer>().material.SetVector("_LightDir", tfLight.forward);
		}
	}
}
