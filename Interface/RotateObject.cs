using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RotateObject : MonoBehaviour {

	public float rotSpeed = 2;

    void OnMouseDrag()
    {
        float rotX = Input.GetAxis("Mouse X")*rotSpeed*Mathf.Deg2Rad;
        //float rotY = Input.GetAxis("Mouse Y")*rotSpeed*Mathf.Deg2Rad;

        transform.RotateAround(Vector3.up, -rotX);
        //transform.RotateAround(Vector3.right, rotY);
    }
}