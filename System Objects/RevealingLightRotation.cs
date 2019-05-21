using UnityEngine;
using System.Collections;


public class RevealingLightRotation : MonoBehaviour {

	Quaternion initRotation;

	void Start () {
		 initRotation = transform.rotation;
	}
	

	void LateUpdate () {
		transform.rotation = initRotation;

	}
}
