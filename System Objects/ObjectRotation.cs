using UnityEngine;
using System.Collections;

public class ObjectRotation : MonoBehaviour {

	public float planetSpeedRotation = 1.0f;

	void Start () {
	
	}
	
	void LateUpdate () {
	
		transform.Rotate(-Vector3.up * Time.deltaTime * planetSpeedRotation);
	}
}
