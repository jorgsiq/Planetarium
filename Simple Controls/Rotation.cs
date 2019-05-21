using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Rotation : MonoBehaviour {

    //velocidade de rotação medida em Xf
    public float speed = 0.5f;

	void Start () {
      
    }
	
	void Update () {

        transform.Rotate(0, speed, 0);
    }
}
