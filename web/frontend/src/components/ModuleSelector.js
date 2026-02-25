import React from 'react';
import { Card, Row, Col, Typography, Badge, List, Button, Space, Tag } from 'antd';
import { CheckCircleOutlined, ArrowRightOutlined } from '@ant-design/icons';

const { Title, Text } = Typography;

function ModuleSelector({ modules, onSelect, uploadedData }) {
  const getStatusColor = (status) => {
    switch (status) {
      case 'ready': return 'success';
      case 'beta': return 'processing';
      case 'dev': return 'default';
      default: return 'default';
    }
  };

  const getStatusText = (status) => {
    switch (status) {
      case 'ready': return '可用';
      case 'beta': return '测试版';
      case 'dev': return '开发中';
      default: return '未知';
    }
  };

  return (
    <div>
      <Title level={4}>🎯 选择分析模块</Title>
      
      {uploadedData && (
        <Card size="small" style={{ marginBottom: 16, background: '#f6ffed' }}>
          <Space>
            <Text>已上传数据:</Text>
            <Tag color="blue">{uploadedData.filename}</Tag>
            <Text>{uploadedData.samples} 样本 × {uploadedData.features} 特征</Text>
          </Space>
        </Card>
      )}

      <Row gutter={[16, 16]}>
        {modules.map(module => (
          <Col xs={24} sm={12} lg={8} key={module.id}>
            <Card
              hoverable
              style={{ height: '100%' }}
              actions={[
                <Button 
                  type="primary" 
                  icon={<ArrowRightOutlined />}
                  onClick={() => onSelect(module.id)}
                >
                  选择此模块
                </Button>
              ]}
            >
              <div style={{ textAlign: 'center', marginBottom: 16 }}>
                <span style={{ fontSize: 48 }}>{module.icon}</span>
                <Title level={5} style={{ marginTop: 8 }}>{module.name}</Title>
              </div>
              
              <List
                size="small"
                dataSource={module.features}
                renderItem={item => (
                  <List.Item>
                    <CheckCircleOutlined style={{ color: '#52c41a', marginRight: 8 }} />
                    {item}
                  </List.Item>
                )}
              />
            </Card>
          </Col>
        ))}
      </Row>
    </div>
  );
}

export default ModuleSelector;
